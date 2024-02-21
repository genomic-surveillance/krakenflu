import pytest
import os.path
import re
from importlib_resources import files

from kraken_flu.src.kraken_db_builder import KrakenDbBuilder

FIXTURE_DIR = files('kraken_flu.tests.fixtures')
TAX_DIR = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','taxonomy'))
FASTA_FILE = FIXTURE_DIR.joinpath(os.path.join('downsampled_ncbi_all_flu','all_influenza_1pcerent_plus_refseq_plus_completeflu.fasta'))

def test_init():
    kdb = KrakenDbBuilder( taxonomy_path= TAX_DIR, fasta_file_path= FASTA_FILE)
    assert kdb
    assert kdb._fasta_handler
    assert kdb._taxonomy_handler
    assert kdb.tax_ids_updated == False
    
def test_update_fasta_tax_ids():
    kdb = KrakenDbBuilder( taxonomy_path= TAX_DIR, fasta_file_path= FASTA_FILE)
    assert kdb.update_fasta_tax_ids(), 'process runs without errors'
    
def test_create_db_ready_dir( tmp_path ):
    kdb = KrakenDbBuilder( taxonomy_path= TAX_DIR, fasta_file_path= FASTA_FILE)

    out_dir = tmp_path / 'out_dir'
    assert not out_dir.is_dir() , 'before we begin, the output dir does not exist'
    kdb.create_db_ready_dir( path= out_dir )
    assert out_dir.is_dir() , 'after calling create_db_ready_dir, the output dir has been created'
    assert os.path.isfile(os.path.join( out_dir, 'taxonomy', 'nodes.dmp')), 'the nodes.dmp file has been created'
    assert os.path.isfile(os.path.join( out_dir, 'taxonomy', 'names.dmp')), 'the names.dmp file has been created'
    assert os.path.isfile(os.path.join( out_dir, 'library', 'library.fna')), 'the FASTA file has been created'
