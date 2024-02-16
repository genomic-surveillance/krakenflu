from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import csv
import shutil
import os.path
import logging

from kraken_flu.src.fasta_handler import FastaHandler
from kraken_flu.src.taxonomy_handler import TaxonomyHandler

logging.basicConfig( format='%(asctime)s %(message)s', level=logging.DEBUG )

class KrakenDbBuilder():   
    """
    This class orchestrates the modifications needed in order to re-organise the influenza 
    genomes in the taxonomy. It uses the FastaHandler and TaxonomyHandler classes for the 
    heavy lifting.
    
    Use this class to create all the files that are needed to build a kraken2 DB with reorganised flu
    genomes.
    
    Consider this example
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
            path to the taxonomy directory from kraken2-build process (kraken2 default name: taxonomy), 
            which must contain the files nodes.dmp and names.dmp
            
        fasta_file_path: str, required
            path to the library (genome data) directory from kraken2-build process (kraken2 default name: library)
        
    """
    
    def __init__( self, taxonomy_path: str, fasta_file_path: str ):
        self.taxonomy_path = taxonomy_path
        self.fasta_file_path = fasta_file_path
    
        self._fasta_handler = FastaHandler( self.fasta_file_path )
        self._taxonomy_handler = TaxonomyHandler( self.taxonomy_path )
        
        self.tax_ids_updated = False
        
    def update_fasta_tax_ids( self ):
        """
        Combines the TaxonomyHandler and FastaHandler to update taxon IDs for output of the new
        FASTA file.
        
        Paramters:
            none
            
        Returns:
            True on success
            
        Side effects:
            - updates taxids in slef._fasta_handler.data
            - sets self.tax_ids_update to True
        """
        if self.tax_ids_updated:
            # already done this, nothing to do
            return True
                
        logging.info( f'starting to assign new taxonomy IDs to sequence records')
        
        # trigger the cascade to create the new taxa
        self._taxonomy_handler.create_influenza_isolate_segment_taxa()
        
        for fasta_record in self._fasta_handler.data:
            # currently, we only make changes to flu sequences
            if fasta_record.is_flu and fasta_record.flu_name and fasta_record.flu_seg_num:
                try:
                    taxid = self._taxonomy_handler.influenza_isolate_segment_tax_ids[ fasta_record.flu_name ][ fasta_record.flu_seg_num ]
                    fasta_record.taxid = taxid
                except KeyError:
                    logging.info( f'found influenza record for isolate "{fasta_record.flu_name}" segment { fasta_record.flu_seg_num} in FASTA file but no corresponding entry in taxonomy - skipping')

        logging.info( f'finished assigning new taxonomy IDs')

        self.tax_ids_updated = True
        return True

    def create_db_ready_dir( self, path: str, force: bool = True, fasta_file_name:str = 'library.fna' ):
        """
        Create a directory with all files needed to build a new kraken database. 
        This method triggers the data acquisition from the FastaHandler and TaxonomyHandler, then runs the
        changing of taxonomy IDs (delegating to the two handlers) before writing a new FASTA file with updated 
        kraken:taxid tags as well as the names and nodes files of the NCBI taxonomy with the added taxa.
        
        The FASTA file is named library.fna as per kraken2 convention but can be given a different name. The 
        NCBI taxonomy files names.dmp and nodes.dmp cannot be renamed as kraken2 relies on those names.
        
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
        
        fasta_file_path = os.path.join( library_path, fasta_file_name )
        names_file_path = os.path.join( taxonomy_path, 'names.dmp')
        nodes_file_path = os.path.join( taxonomy_path, 'nodes.dmp')

        if os.path.exists( path ): 
            if not force:
                raise ValueError(f'directory { path } exists already. Will not write into existing directory')
        else:
            os.mkdir( path )
            
        if not os.path.exists( library_path ):
            os.mkdir( library_path )
        if not os.path.exists( taxonomy_path ):
            os.mkdir( taxonomy_path )
        
        logging.info( f'found { self._fasta_handler.n_seq_total} sequence records in {self.fasta_file_path}')
        logging.info( f'{ self._fasta_handler.n_seq_flu} sequence records identified as influenza')
        self._fasta_handler.remove_incomplete_flu()
        logging.info( f'{ self._fasta_handler.n_seq_filtered} sequence records have been removed as incomplete flu genomes')
        
        if not self.tax_ids_updated:
            self.update_fasta_tax_ids()
        
        logging.info( f'writing output files')
        self._fasta_handler.write_fasta( fasta_file_path )
        self._taxonomy_handler.write_names_file( names_file_path )
        self._taxonomy_handler.write_nodes_file( nodes_file_path )
        logging.info( f'---- process complete ----')
        return True

