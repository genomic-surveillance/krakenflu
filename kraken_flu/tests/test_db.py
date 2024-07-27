import pytest
import os.path
from importlib_resources import files
from sqlalchemy import select

from kraken_flu.src.db import Db, Sequence

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
    rows = db._session.execute(select(Sequence.dna_sequence).where( Sequence.fasta_header=='some test fasta header')).all()
    assert len(rows) == 0, 'before inserting, we have no result'
    
    db.add_sequence(
        fasta_header= 'some test fasta header',
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
    # retrieve the newly inserted sequence row
    rows = db._session.execute(select(Sequence.dna_sequence).where( Sequence.fasta_header=='some test fasta header')).all()
    assert len(rows) == 1, 'after inserting, we have a single result'
    assert rows[0].dna_sequence == 'ACGCATCGAA' ,'... the result has the correct DNA sequence value'