from Bio import SeqIO
import re
import csv
import shutil
import os.path
from cached_property import cached_property
import logging

logging.basicConfig( format='%(asctime)s %(message)s', level=logging.DEBUG )

class KrakenDbNcbiFiles():   
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
    KRAKEN_TAX_ID_ASSIGNMENT_REGEX = re.compile(r'kraken:taxid\|([0-9]+)\|')
    NCBI_ACC_REGEX = re.compile(r'[A-Z]{2}_[0-9]{6,}\.[0-9]')
    
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
        regexes = [ 
            re.compile(r'Influenza A.*H1N1'), 
            re.compile(r'Influenza A.*H3N2'),
            re.compile(r'Influenza B') ]
        
        new_tax_id_i = self.min_new_tax_id
        data = {}
        need_acc2tax_scan = False
        logging.info( f'scanning file { self.fasta_file_path } for influenza viruses')
        flu_count = 0
        with open( self.fasta_file_path ) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                for regex in regexes:
                    if regex.search( record.description ):
                        flu_count += 1
                        # the currently assigned tax ID becomes the new parent tax ID
                        parent_tax_id = self._parse_kraken_tax_id( record.description )
                        if parent_tax_id is None:
                            need_acc2tax_scan = True
                        
                        ncbi_id = self._parse_ncbi_accession_id( record.description )
                        if ncbi_id in data:
                            raise ValueError(f'NCBI ID {ncbi_id} found more than once in FASTA file {self.fasta_file_path}' )
                        
                        # clean up the FASTA header to serve as a name for this segment
                        mod_name = self.KRAKEN_TAX_ID_ASSIGNMENT_REGEX.sub( '', record.description )
                        
                        data[ ncbi_id ] = {
                            'name': mod_name, 
                            'new_tax_id': new_tax_id_i,
                            'new_parent_id': parent_tax_id }
                                                
                        new_tax_id_i += 1
                        continue
        
        logging.info( f'done - found { flu_count } segment sequences in FASTA file')

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
                ( _, acc_version, taxid, _ ) = row
                if acc_version in data:
                    if data[ acc_version ]['new_parent_id'] is None:
                        data[ acc_version ]['new_parent_id'] = taxid
                        
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
        
        try: 
            with open( path, 'w' ) as out_fh:
                try:
                    with open( self.fasta_file_path ) as in_fh:
                        logging.info(f'writing modified FASTA file to { path }')
                        for record in SeqIO.parse( in_fh, "fasta"):
                            ncbi_id = self._parse_ncbi_accession_id( record.description )
                            existing_kraken_tax_id = self._parse_kraken_tax_id( record.description )
                            if ncbi_id in self.flu_genomes_ncbi_to_new_tax_and_parent_ids:
                                # this is a segmented flu genome
                                new_kraken_tax_id = self.flu_genomes_ncbi_to_new_tax_and_parent_ids[ ncbi_id ]['new_tax_id']
                                new_kraken_tax_str = f'kraken:taxid|{new_kraken_tax_id}|'
                                new_desc = self.KRAKEN_TAX_ID_ASSIGNMENT_REGEX.sub( new_kraken_tax_str, record.description )
                                record.description = new_desc
                            SeqIO.write( record, out_fh, "fasta" )   
                except (PermissionError, FileNotFoundError) as e:
                    raise ValueError(f'cannot read from { self.fasta_file_path }')
                
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
            for ncbi_id in self.flu_genomes_ncbi_to_new_tax_and_parent_ids.keys():
                name = self.flu_genomes_ncbi_to_new_tax_and_parent_ids[ ncbi_id ]['name']
                tax_id = self.flu_genomes_ncbi_to_new_tax_and_parent_ids[ ncbi_id ]['new_tax_id']
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
            for ncbi_id in self.flu_genomes_ncbi_to_new_tax_and_parent_ids.keys():
                tax_id = self.flu_genomes_ncbi_to_new_tax_and_parent_ids[ ncbi_id ]['new_tax_id']
                parent_tax_id = self.flu_genomes_ncbi_to_new_tax_and_parent_ids[ ncbi_id ]['new_parent_id']

                # the hardcoded '9' in field 5 is for virus division of the NCBI taxonomy
                out_fh.write(
                    "\t|\t".join( [ str( tax_id ), str( parent_tax_id ), 'no rank', '', '9', '1', '1', '1', '0', '1', '1', '', '' ]) + "\t|\n")

        logging.info('finished writing modified nodes file')
        return True
    
    def write_prelim_map_file( self, path: str ):
        """
        Writes the prelim_map.txt file (if one exists) with the new modified taxon IDs for 
        the segmented virus but leaving other IDs unchanged.
                
        Parameters:
            path: str
                output file full path
                
        Returns:
            True if success
            
        Side effects:
            Writes to a file 'path'
            
        """
        if not self.prelim_map_file_path:
            logging.info(f'no prelim_map.txt file provided in inputs - nothing to do')
            return True
        
        try: 
            with open( path, 'w' ) as out_fh:
                try:
                    with open( self.prelim_map_file_path ) as in_fh:
                        logging.info(f'writing modified prelim_map.txt file to { path }')
                        for row in in_fh:
                            row = row.rstrip()
                            fields = row.split("\t")
                            ncbi_id = self._parse_ncbi_accession_id( fields[1] )
                            if ncbi_id in self.flu_genomes_ncbi_to_new_tax_and_parent_ids:
                                # this is a segmented flu genome
                                new_kraken_tax_id = self.flu_genomes_ncbi_to_new_tax_and_parent_ids[ ncbi_id ]['new_tax_id']
                                fields[1] = f'kraken:taxid|{new_kraken_tax_id}|{ncbi_id}'
                                fields[2] = str(new_kraken_tax_id)
                            
                            out_fh.write( "\t".join( fields ) + "\n" )
                except (PermissionError, FileNotFoundError) as e:
                    raise ValueError(f'cannot read from { self.fasta_file_path }')
                
        except (PermissionError, FileNotFoundError) as e:
            raise ValueError(f'path { path } does not exist or is not writable, cannot create file')

        logging.info('finished writing modified prelim_map.txt file')
        return True
    
    def create_db_ready_dir( self, path: str ):
        """
        Create a directory with all files needed to build a new kraken database. This includes the files that
        are copied unchanged as well as the files that are modified, which are:
            - fasta file
            - prelim_map.txt (preliminary mapping of NCBI ID to taxon ID, not always present)
            - names.dmp
            - nodes.dmp
            
        If a NCBI acc2taxid file was provided, it is also copied across. Beware that this is a large file (currently around
        12GB) even if you are building a small database because it contains an entry for every record in NCBI.
        
        The final directory can be used by kraken2 build as follows:
        
        kraken2-build --build --db PATH
        
        Parameters:
            path: str
                path to the directory into which the files will be written.
                This directory must NOT exist yet. This is a precaution to ensure that the user does
                not mix files for DB building which can result in an inconsistent built.
                
        Returns:
            True if success
            
        Side effects:
            Writes files to path
            
        """
        if os.path.exists( path ):
            raise ValueError(f'directory { path } exists already. Will not write into existing directory')
        os.mkdir( path )
        
        library_path = os.path.join( path , 'library' )
        taxonomy_path = os.path.join( path, 'taxonomy' )

        os.mkdir( library_path )
        os.mkdir( taxonomy_path )
    
        self.write_modified_fasta_file( os.path.join( library_path, 'library.fna' ))
        self.write_prelim_map_file( os.path.join( library_path, 'prelim_map.txt' ))
        self.write_modified_names_files( os.path.join( taxonomy_path, 'names.dmp' ))
        self.write_modified_nodes_files( os.path.join( taxonomy_path, 'nodes.dmp' ))

        if self.acc2tax_file_path:
            shutil.copyfile( self.nodes_file_path, os.path.join( path, 'taxonomy', 'nucl_gb.accession2taxid' ) )
    
    def _parse_ncbi_accession_id( self, in_str: str):
        """
        Parse the NCBI accession ID from a string (like FASTA header)
        
        Parameters:
            in_str: str, required
                the string that will be parsed
                
        Returns:
            NCBI accession ID (int)
            Throws exception if not present
        """
        try:
            return self.NCBI_ACC_REGEX.search( in_str ).group(0)
        except AttributeError:
            raise ValueError(f"could not parse NCBI acc ID from '{in_str}'")

    def _parse_kraken_tax_id( self, in_str: str):
        """
        Attempts to parse a taxon ID from the input string, which is present in pre-processed 
        kraken2 database build files. If it is not found, None is returned.
                
        Parameters:
            in_str: str, required
                the string that will be parsed
                
        Returns:
            taxon ID if present (int), otherwise None
        """
        try:
            return int( self.KRAKEN_TAX_ID_ASSIGNMENT_REGEX.search( in_str ).group( 1 ) ) 
        except AttributeError:
            return None