import re
import os.path
from cached_property import cached_property
import logging

logging.basicConfig( format='%(asctime)s %(message)s', level=logging.DEBUG )

class TaxonomyHandler():
    """
    This class handles the NCBI taxonomy, its main task is to restructure the taxonomy for
    influenza isolates.
    
    Parameters:
        taxonomy_path: str, required
            Path to the taxonomy directory from kraken2-build process (kraken2 default name: taxonomy)
            It needs to contain the names.dmp and nodes.dmp files from NCBI
    """
    
    def __init__( self, taxonomy_path: str ):
        
        if not os.path.isdir( taxonomy_path ):
            raise ValueError(f'path { taxonomy_path} does not exist or is not a directory')
        else:
            self.taxonomy_path = taxonomy_path        
        
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
        
        
