import pytest
import os.path
from importlib_resources import files

from kraken_flu.src.taxonomy_loader import load_taxonomy, _load_names, _load_nodes, _read_tax_data_file_row, load_seq2taxid, _get_num_records
from kraken_flu.src.db import Db

FIXTURE_DIR = files('kraken_flu.tests.fixtures')
NODES_FILE = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','taxonomy','nodes.dmp'))
NAMES_FILE = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','taxonomy','names.dmp'))
ACC2TAXID_FILE = FIXTURE_DIR.joinpath(os.path.join('ncbi_data_not_kraken','nucl_gb.accession2taxid'))

DEBUG=True
#DEBUG=False

@pytest.fixture(scope='function')
def setup_db( tmp_path ):
    """ initialise the DB """
    db_path = tmp_path / 'kraken_flu.db'
    db = Db(db_path, debug=DEBUG)
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
    rows = db._cur.execute("select * FROM taxonomy_names").fetchall()
    assert len(rows) == 0, 'no taxonomy_names data in DB before we start uploading'
    assert _load_names(db=db, names_file_path=NAMES_FILE), '_load_names returns True'
    
    rows = db._cur.execute("select * FROM taxonomy_names").fetchall()
    assert len(rows) == 52, 'after running _load_names, we now have 20 rows of data in taxonomy_names table'

def test__load_nodes(setup_db):
    db = setup_db
    # NOTE: the db session should not normally be accessed directly. Only doing this here so we can perform a test 
    # that doesn't rely on a Db class method. 
    rows = db._cur.execute("select * FROM taxonomy_nodes").fetchall()
    assert len(rows) == 0, 'no taxonomy_nodes data in DB before we start uploading'
    assert _load_nodes(db=db, nodes_file_path=NODES_FILE), '_load_nodes returns True'
    
    rows = db._cur.execute("select * FROM taxonomy_nodes").fetchall()
    assert len(rows) == 20, 'after running _load_names, we now have 20 rows of data in taxonomy_nodes table'
    
def test_load_taxonomy(setup_db):
    db=setup_db
    assert len(db._cur.execute("select * FROM taxonomy_nodes").fetchall()) ==  len(db._cur.execute("select * FROM taxonomy_names").fetchall()) == 0, 'before we start uploading, no taxonomy nodes or names are in the DB'
    assert load_taxonomy(db=db, names_file_path=NAMES_FILE, nodes_file_path=NODES_FILE)
    assert len(db._cur.execute("select * FROM taxonomy_nodes").fetchall()) == 20, 'after running load_taxonomy 20 taxonomy_nodes exist in the DB'
    assert len(db._cur.execute("select * FROM taxonomy_names").fetchall()) == 52, 'after running load_taxonomy 52 taxonomy_names exist in the DB'

    # retrieve a node and corresponding names by tax_id and name class 
    rows = db._cur.execute(
        """
        SELECT taxonomy_nodes.tax_id, taxonomy_nodes.parent_tax_id
        FROM taxonomy_nodes 
        INNER JOIN taxonomy_names 
        WHERE taxonomy_nodes.tax_id == ? AND taxonomy_names.name_class == ?
        GROUP BY taxonomy_nodes.tax_id, taxonomy_nodes.parent_tax_id
        """,
        (2697049, 'scientific name',
        )).fetchall()
    assert len(rows) == 1, 'a single result is returned'
    assert rows[0]['parent_tax_id'] == 694009, 'parent tax_id is correct'
    #assert [ r.name for r in row.TaxonomyNode.taxonomy_names ] == ['Severe acute respiratory syndrome coronavirus 2'], 'a single name is retrieved when filtering for scientific name'
    
    def test__load_seq2taxids(setup_db):
        db = setup_db
        assert len(db._cur.execute("select * FROM seq2taxid").fetchall()) ==  len(db._cur.execute("select * FROM seq2taxid").fetchall()) == 0, 'before we start uploading, no seq2taxid records are present in the DB'
        assert load_seq2taxids(db=db, acc2taxid_file_path=ACC2TAXID_FILE)
        assert len(db._cur.execute("select * FROM seq2taxid").fetchall()) == 35, 'after running load_seq2taxids 20 accessions exist in the DB'

        #TODO complete test
        
def test__get_num_records():
    assert _get_num_records(NAMES_FILE) == 52, 'there are 52 lines in the names file'
    
