from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import csv
import shutil
import os.path
from cached_property import cached_property
import logging

from kraken_flu.src.fasta_parser import FastaParser

logging.basicConfig( format='%(asctime)s %(message)s', level=logging.DEBUG )

class KrakenDbBuilder():   
    """
    This class handles the file collection that is created from the kraken2 build process.
    The class provides methods for the replacement of all influenza whole genomes with individual
    segmented references.
    
    To illustrate the meaning of the required attributes, consider this example
    of a library build on NCBI viral refseq:
    
    The kraken2 DB preparation commands are:
    kraken2-build \
        --download-taxonomy \
        --db SOME/PATH/
    
    kraken2-build \
        --download-library viral \
        --db SOME/PATH/
        
    This creates the following files in SOME/PATH
    ├── library
    │   └── viral
    │       ├── assembly_summary.txt
    │       ├── library.fna
    │       ├── library.fna.masked
    │       ├── manifest.txt
    │       └── prelim_map.txt
    └── taxonomy
        ├── accmap.dlflag
        ├── citations.dmp
        ├── delnodes.dmp
        ├── division.dmp
        ├── gc.prt
        ├── gencode.dmp
        ├── images.dmp
        ├── merged.dmp
        ├── names.dmp
        ├── nodes.dmp
        ├── nucl_gb.accession2taxid
        ├── nucl_wgs.accession2taxid
        ├── readme.txt
        ├── taxdump.dlflag
        ├── taxdump.tar.gz
        └── taxdump.untarflag
        
    The two required attributes are the full path to 'library' and 'taxonomy'. 
    The files that are read and manipulated are:
    library/library.fna: sequence data
    library/prelim_map.txt: mapping of sequence IDs to taxon IDs
    taxonomy/names.dmp: lists all taxonomy names, will have new names added for flu segment
    taxonomy/nodes.dmp: lists all nodes in the taxonomy, will have new nodes added for flu segments.
    
    When not using pre-built kraken2 reference files, a mapping from NCBI accession ID to taxonomy ID is
    required. This can be obtained from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
    
    Parameters:
        taxonomy_path: str, required
            path to the taxonomy directory from kraken2-build process (kraken2 default name: taxonomy)
            
        library_path: str, required
            path to the library (genome data) directory from kraken2-build process (kraken2 default name: library)
            
        acc2tax_file_path: str, optional
            path to the NCBI file nucl_gb.accession2taxid if one is used (only relevant when not using kraken2 
            pre-built DB reference files)
        
    """
    
    def __init__( self, taxonomy_path: str, library_path: str, acc2tax_file_path: str = None ):
        self.taxonomy_path = taxonomy_path
        self.library_path = library_path
        self.acc2tax_file_path = acc2tax_file_path
        
        if not os.path.isdir( self.taxonomy_path ):
            raise ValueError(f'path { self.taxonomy_path} does not exist or is not a directory')
        if not os.path.isdir( self.library_path ):
            raise ValueError(f'path { self.library_path} does not exist or is not a directory')
        
        names_file_path = os.path.join(self.taxonomy_path, 'names.dmp')
        if os.path.exists( names_file_path ):
            self.names_file_path = names_file_path
        else:
            raise ValueError(f'missing file { names_file_path }')

        nodes_file_path = os.path.join(self.taxonomy_path, 'nodes.dmp')
        if os.path.exists( nodes_file_path ):
            self.nodes_file_path = nodes_file_path
        else:
            raise ValueError(f'missing file { nodes_file_path }')

        fasta_file_path = os.path.join(self.library_path, 'library.fna')
        if os.path.exists( fasta_file_path ):
            self.fasta_file_path = fasta_file_path
        else:
            raise ValueError(f'missing file { fasta_file_path }')
        
        if self.acc2tax_file_path and not os.path.isfile( self.acc2tax_file_path ):
            raise ValueError(f'file { self.acc2tax_file_path } does not exists')
        
        prelim_path =  os.path.join( self.library_path, 'prelim_map.txt')
        if os.path.exists( prelim_path ):
            self.prelim_map_file_path = prelim_path
        else:
            self.prelim_map_file_path = None
        
        self.min_new_tax_id = self._find_max_used_tax_id() + 1
        
    def _find_max_used_tax_id( self ):
        """
        Returns the maximum taxonomy ID found in the names file
        """
        f = open( self.names_file_path, 'r' )
        max_id = 0
        for line in f:
            id = int(line.split("\t|\t")[0])
            if id > max_id:
                max_id = id
                
        return max_id 
    
    @cached_property
    def fasta_data(self):
        """
        The FASTA header and seq data, delegated to FastaParser
        """
        fp = FastaParser( self.fasta_file_path )
        return fp.data
    
    @cached_property
    def flu_genomes_ncbi_to_new_tax_and_parent_ids( self ):
        """
        Builds a dictionary of NCBI IDs to new (created by this utility) taxonomy ID and a parent
        tax ID which is the old taxonomy ID of this record.
        The information is obtained by reading the FASTA file library.fna and, if taxonomy Ids are not in that
        file (not a prebuilt KRAKEN2 DB file) the NCBI file nucl_gb.accession2taxid is used to obtain the mapping
        between NCBI accession IDs and taxonomy IDs
        
        Example:
        
        FASTA headers of two of the segments of a given isolate:
        >kraken:taxid|211044|NC_002023.1 Influenza A virus (A/Puerto Rico/8/1934(H1N1)) segment 1, complete sequence
        >kraken:taxid|211044|NC_002021.1 Influenza A virus (A/Puerto Rico/8/1934(H1N1)) segment 2, complete sequence

        These are from a KRAKEN2 prebuilt DB file, so they already contain the tax ID in the header as
        kraken:taxid
        Both have the same tax_id 211044.
        
        This methods assigns each NCBI ID a new "made up" tax ID and both are assigned tax ID 211044 as the parent 
        taxonomy ID. The parent ID will later be recorded in the modified nodes.dmp file. As a result, all segments are still
        grouped as one genome under the tax ID 211044 but reads can be assigned to a particular segment in the final report of KRAKEN2. 
        
        FASTA header description and the segment number.
        Reads the FASTA file an retrieves all relevant flu genome headers.
        For each NCBI ID found in a flu sequence header, a new Tax Id is created and stored in
        the resulting dictionary. 
        This does not guarantee that segments of a single genome get assigned subsequent tax ID 
        numbers but that does not matter for the database we are building.
        
        Parameters:
            none
            
        Returns:
            dict with the following structure:
                { NCBI accession ID 1: 
                    {   new_tax_id: new tax ID1, 
                        new_parent_id: parent tax ID 1 }, ... }
        
        """
        
        new_tax_id_i = self.min_new_tax_id
        data = {}
        need_acc2tax_scan = False
        logging.info( f'scanning file { self.fasta_file_path } for influenza viruses')
        
        fasta_data = self.fasta_data
        for record in fasta_data:
            if record.is_flu:
                parent_tax_id = record.taxid
                if parent_tax_id is None:
                    need_acc2tax_scan = True
                
                ncbi_id = record.ncbi_acc
                data[ ncbi_id ] = {
                    'new_tax_id': new_tax_id_i,
                    'new_parent_id': parent_tax_id }
                                        
                new_tax_id_i += 1
                continue

        if need_acc2tax_scan:
            logging.info('reading NCBI acc2taxid file to assign taxon IDs to NCBI IDs')
            if not self.acc2tax_file_path:
                # TODO: could check for a prelim_maps.txt file in the taxonomy dir and use that?
                #       Not sure if there is any merit because if that file exists, it seems the data is
                #       downloaded from kraken2 pre-built files and will have tax IDs in FASTA headers anyway(?)
                raise ValueError('at least some sequences have no tax ID in FASTA header and no NCBI acc2tax file was provided')
            
            data = self._get_tax_ids_from_acc2taxid_map( data )
            logging.info('done reading NCBI acc2taxid file')
        return data
    
    def _get_tax_ids_from_acc2taxid_map( self, data ):
        """
        Reads the large NCBI mapping file from NCBI accession to taxonomy ID once and extracts the 
        missing taxonomy IDs for the data dictionary.

        Parameters:
            data: dict
                data structure keyed by NCBI accession ID
            
        Returns:
            dict of dicts: same as input but with 'new_parent_id' filled in 
            
        """
        
        with open( self.acc2tax_file_path, "r", encoding="utf8") as fh:
            tsv_reader = csv.reader( fh, delimiter="\t")
            # Skip the first row, which is the header
            next(tsv_reader)
            for row in tsv_reader:
                ( acc_id, acc_version, taxid, _ ) = row
                
                # lookup on accession ID or accession ID version (.1, .2, etc)
                # if one of them is present and matches the recorded NCBI IDs, this becomes the
                # new parent ID
                if acc_version in data:
                    use_id = acc_version
                elif acc_id in data:
                    use_id = acc_id
                else:
                    use_id = None
                
                if use_id:
                    if data[ use_id ]['new_parent_id'] is None:
                        data[ use_id ]['new_parent_id'] = int(taxid)
                        
        return data
    
    def write_modified_fasta_file( self, path: str ):
        """
        Writes the FASTA data from the input file to a new file with modified (or newly assigned) kraken2
        taxon IDs in the FASTA headers.
        If the input FASTA file already contains kraken2 tax IDs, these will be overwritten with the newly
        assigned IDs where they exist (ie the segmented flu genomes). 
        Flu genomes that do not yet have a kraken2 taxon ID are now written with this ID in the kraken2 format
        as in this example:
        
        >kraken:taxid|211044|NC_002023.1 Influenza A virus (A/Puerto Rico/8/1934(H1N1))
        
        Parameters:
            path: str
                output file full path
                
        Returns:
            True if success
            
        Side effects:
            Writes to a file 'path'
            
        """
        fasta_data = self.fasta_data
        acc2tax_data = self.flu_genomes_ncbi_to_new_tax_and_parent_ids
        
        try: 
            with open( path, 'w' ) as out_fh:
                logging.info(f'writing modified FASTA file to { path }')
                for record in fasta_data:
                    final_head = record.mod_head
                    if record.is_flu:
                        if record.ncbi_acc in acc2tax_data:
                            new_kraken_tax_id = acc2tax_data[record.ncbi_acc]['new_tax_id']
                            new_kraken_tax_str =  f'kraken:taxid|{new_kraken_tax_id}|'
                            final_head = new_kraken_tax_str + final_head
                        else:
                            raise ValueError(f'record {record.ncbi_acc} is marked as "flu" but has no data in accession2taxid map')
                            
                    sr = SeqRecord(
                        Seq(record.sequence),
                        id= final_head,
                        description= ''
                    )
                    SeqIO.write( sr, out_fh, "fasta" )
        except (PermissionError, FileNotFoundError) as e:
            raise ValueError(f'path { path } does not exist or is not writable, cannot create file')

        logging.info('finished writing modified FASTA file')
        return True
    
    def write_modified_names_files( self, path: str ):
        """
        Writes the names.dmp file with added new taxa from the segmented flu genomes.
                
        Parameters:
            path: str
                output file full path
                
        Returns:
            True if success
            
        Side effects:
            Writes to a file 'path'
            
        """
        try:
            # most of the names file is just copied from the original
            shutil.copyfile( self.names_file_path, path )
        except (PermissionError, FileNotFoundError) as e:
            raise  ValueError(f'could not copy from { self.names_file_path } to { path }: { e }')
        
        # now append the new segmented genome names
        logging.info(f'writing modified names file to { path }')
        with open( path, 'a' ) as out_fh:
            for record in self.fasta_data:
                if record.is_flu and record.ncbi_acc in self.flu_genomes_ncbi_to_new_tax_and_parent_ids:
                    try:
                        name = 'Influenza ' + record.flu_name + ' segment ' + str(record.flu_seg_num)
                    except TypeError:
                        raise ValueError(f'no flu isolate name or segment number recorded for {record.ncbi_acc}')
                    tax_id = self.flu_genomes_ncbi_to_new_tax_and_parent_ids[record.ncbi_acc]['new_tax_id']
                    if not tax_id:
                        raise ValueError(f'no taxid recorded for {record.ncbi_acc}')
                    out_fh.write( 
                        "\t|\t".join( [ str( tax_id ) , str( name ), '', 'scientific name'  ] )  + "\t|\n" )

        logging.info('finished writing modified names file')
        return True
    
    def write_modified_nodes_files( self, path: str ):
        """
        Writes the nodes.dmp file with added new taxa from the segmented flu genomes.
                
        Parameters:
            path: str
                output file full path
                
        Returns:
            True if success
            
        Side effects:
            Writes to a file 'path'
            
        """
        try:
            # most of the names file is just copied from the original
            shutil.copyfile( self.nodes_file_path, path )
        except (PermissionError, FileNotFoundError) as e:
            raise  ValueError(f'could not copy from { self.nodes_file_path } to { path }: { e }')
        
        # now append the new segmented genome names
        logging.info(f'writing modified nodes file to { path }')
        with open( path, 'a' ) as out_fh:
            for record in self.fasta_data:
                if record.is_flu and record.ncbi_acc in self.flu_genomes_ncbi_to_new_tax_and_parent_ids:
                    tax_id = self.flu_genomes_ncbi_to_new_tax_and_parent_ids[record.ncbi_acc]['new_tax_id']
                    parent_tax_id = self.flu_genomes_ncbi_to_new_tax_and_parent_ids[record.ncbi_acc]['new_parent_id']
                    if not tax_id:
                        raise ValueError(f'no taxid recorded for {record.ncbi_acc}')
                    if not parent_tax_id:
                        raise ValueError(f'no parent taxid recorded for {record.ncbi_acc}')
                    
                    # the hardcoded '9' in field 5 is for virus division of the NCBI taxonomy
                    out_fh.write(
                        "\t|\t".join( [ str( tax_id ), str( parent_tax_id ), 'no rank', '', '9', '1', '1', '1', '0', '1', '1', '', '' ]) + "\t|\n")

        logging.info('finished writing modified nodes file')
        return True
    
    def create_db_ready_dir( self, path: str, force: bool = True ):
        """
        Create a directory with all files needed to build a new kraken database. This includes the files that
        are copied unchanged as well as the files that are modified, which are:
            - fasta file (will ba saved as library.fna to comply with kraken2 convention)
            - names.dmp
            - nodes.dmp
            
        The final directory can be used by kraken2 build as follows.
        Assuming the db-ready dir (path parameter) is /path/to/dbready_dir/
        and the path to the DB we are building is /path/to/new_db/
        
        1) add to DB taxonomy resources
        this will overwrite the exiting names and nodes file in the DB taxonomy folder, which is ok
        because the modified files contain only additional content
        
        $ cp /path/to/dbready_dir/*dmp /path/to/new_db/taxonomy
        
        2) add to sequence library
        use kraken2-build to ad the sequence FASTA file to the DB resources
        
        $ kraken2-build \
            --add-to-library /path/to/dbready_dir/library.fna \
            --db  /path/to/new_db/
        
        3) build the DB
        kraken2-build --build --db  /path/to/new_db/
        
        
        Parameters:
            path: str
                path to the directory into which the files will be written.
                If force is not used, the DIR must not exist yet
                
            force: bool, optional, defaults to True
                if true, force overwrite whatever is already in the output dir. If not, will fail
                if directory already exists
                
        Returns:
            True if success
            
        Side effects:
            Writes files to path
            
        """
        library_path = os.path.join( path , 'library' )
        taxonomy_path = os.path.join( path, 'taxonomy' )
        
        if os.path.exists( path ): 
            if not force:
                raise ValueError(f'directory { path } exists already. Will not write into existing directory')
        else:
            os.mkdir( path )
            
        if not os.path.exists( library_path ):
            os.mkdir( library_path )
        if not os.path.exists( taxonomy_path ):
            os.mkdir( taxonomy_path )
    
        self.write_modified_fasta_file( os.path.join( library_path, 'library.fna' ))
        self.write_modified_names_files( os.path.join( taxonomy_path, 'names.dmp' ))
        self.write_modified_nodes_files( os.path.join( taxonomy_path, 'nodes.dmp' ))
