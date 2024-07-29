import pytest
import os.path
from importlib_resources import files
from sqlalchemy import select

from kraken_flu.src.taxonomy_loader import load_taxonomy, _load_names, _load_nodes, _read_tax_data_file_row
from kraken_flu.src.db import Db, TaxonomyName, TaxonomyNode

FIXTURE_DIR = files('kraken_flu.tests.fixtures')
NODES_FILE = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','taxonomy','nodes.dmp'))
NAMES_FILE = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','taxonomy','names.dmp'))

@pytest.fixture(scope='function')
def setup_db( tmp_path ):
    """ initialise the DB """
    db_path = tmp_path / 'kraken_flu.db'
    db = Db(db_path)
    yield db
    
def test__read_tax_data_file_row():
    row='11320	|	2955291	|	no rank	|		|	9	|	1	|	1	|	1	|	0	|	1	|	0	|	0	|	code compliant; specified	|'
    expected_result = [
        '11320','2955291','no rank','','9','1','1','1','0','1','0','0','code compliant; specified']
    d=_read_tax_data_file_row(row)
    assert d == expected_result, 'row parsed correctly'

    
def test__load_names(setup_db):
    db = setup_db
    # NOTE: the db session should not normally be accessed directly. Only doing this here so we can perform a test 
    # that doesn't rely on a Db class method. 
    rows = db._session.execute(select(TaxonomyName)).all()
    assert len(rows) == 0, 'no taxonomy_names data in DB before we start uploading'
    assert _load_names(db=db, names_file_path=NAMES_FILE), '_load_names returns True'
    
    rows = db._session.execute(select(TaxonomyName)).all()
    assert len(rows) == 52, 'after running _load_names, we now have 20 rows of data in taxonomy_names table'

def test__load_nodes(setup_db):
    db = setup_db
    # NOTE: the db session should not normally be accessed directly. Only doing this here so we can perform a test 
    # that doesn't rely on a Db class method. 
    rows = db._session.execute(select(TaxonomyNode)).all()
    assert len(rows) == 0, 'no taxonomy_nodes data in DB before we start uploading'
    assert _load_nodes(db=db, nodes_file_path=NODES_FILE), '_load_nodes returns True'
    
    rows = db._session.execute(select(TaxonomyNode)).all()
    assert len(rows) == 20, 'after running _load_names, we now have 20 rows of data in taxonomy_nodes table'