import pytest
import os.path
import re
from importlib_resources import files

from kraken_flu.src.taxonomy_handler import TaxonomyHandler

FIXTURE_DIR = files('kraken_flu.tests.fixtures')
TAX_DIR = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','taxonomy'))

def test_init():
    th = TaxonomyHandler( taxonomy_path= TAX_DIR)
    assert th
    
    assert os.path.basename( th.names_file_path ) == 'names.dmp', 'found the names file'
    assert os.path.basename( th.nodes_file_path ) == 'nodes.dmp', 'found the nodes file'
    
    
