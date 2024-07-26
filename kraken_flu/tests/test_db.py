import pytest
import os.path
from importlib_resources import files

from kraken_flu.src.db import Db

@pytest.fixture(scope='function')
def setup_db( tmp_path ):
    db_path = tmp_path / 'kraken_flu.db'
    db = Db(db_path)
    yield db
    
def test_db_connect( tmp_path ):
    db_path = tmp_path / 'kraken_flu.db'
    assert not os.path.exists(db_path), 'before we build a DB, the DB file does not exist'
    db = Db(db_path)
    assert os.path.exists(db_path), 'after creating a DB, the DB file exists now'

def test_add_sequence(setup_db):
    db = setup_db
    db.add_sequence(
        fasta_header= 'some fasta header',
        dna_sequence= 'ACGCATCGAA',
        category= None,
        flu_type= 'A',
        ncbi_acc= None,
        original_taxid= 12345,
        is_flu= True,
        isolate_name= 'A/name/of/flu',
        segment_number= 2, 
        h_subtype= 3,
        n_subtype= 2
    )
    