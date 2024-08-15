import pytest
import os.path
from importlib_resources import files

from kraken_flu.src.db import Db

# set this to True to make the tests print all executed SQL or False to stop that
PRINT_SQL_TRACE=False

@pytest.fixture(scope='function')
def setup_db( tmp_path ):
    db_path = tmp_path / 'kraken_flu.db'
    db = Db(db_path, debug=PRINT_SQL_TRACE)
    yield db

@pytest.fixture(scope='function')
def setup_db_with_fixture( setup_db ):
    db = setup_db
    fixture_dir = files('kraken_flu.tests.fixtures')
    db_fixture_file = fixture_dir.joinpath(os.path.join('db_fixtures','db_fixture1.sql'))
    with open(db_fixture_file, 'r') as file:
        file_content = file.read()
    db._cur.executescript(file_content)
    yield db
    
@pytest.fixture(scope='function')
def setup_db_with_real_world_fixture( setup_db ):
    db = setup_db
    fixture_dir = files('kraken_flu.tests.fixtures')
    db_fixture_file = fixture_dir.joinpath(os.path.join('db_fixtures','db_fixtures2.sql'))
    with open(db_fixture_file, 'r') as file:
        file_content = file.read()
    db._cur.executescript(file_content)
    yield db
    
def test_db_connect( tmp_path ):
    db_path = tmp_path / 'kraken_flu.db'
    assert not os.path.exists(db_path), 'before we build a DB, the DB file does not exist'
    db = Db(db_path)
    assert os.path.exists(db_path), 'after creating a DB, the DB file exists now'

def test_add_sequence(setup_db):
    db = setup_db
    stmt = "SELECT * FROM sequences WHERE fasta_header == 'some test fasta header'"
    rows = db._cur.execute(stmt).fetchall()
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
    rows = db._cur.execute(stmt).fetchall()
    assert len(rows) == 1, 'after inserting, we have a single result'
    row = rows[0]
    assert row['dna_sequence'] == 'ACGCATCGAA' ,'... the result has the correct DNA sequence value'
    
    # We record the original kraken:taxid if there is one in the header but this does not assign
    # the record to a node in the taxonomy_nodes table and hence does not assign the tax_id field
    # Instead, we record any kraken:taxid from the FASTA header in the original_tax_id field
    assert row['tax_id'] == None, '...no taxonomy node has been assigned yet, so not tax_id is assigned'
    assert row['original_tax_id'] == 12345, '... the original kraken taxid from the inserted record is correct'

def test_get_leaf_node_tax_ids(setup_db_with_fixture):
    db = setup_db_with_fixture
    leaf_node_tax_ids = db.get_leaf_node_tax_ids()
    assert isinstance( leaf_node_tax_ids, list)
    assert len(leaf_node_tax_ids) == 2, 'the fixtures have 2 leaf nodes'
    assert set(leaf_node_tax_ids) == set([3,4]), 'correct tax_ids identified as leaf nodes'
    
def test_get_parent_tax_id(setup_db_with_fixture):
    db = setup_db_with_fixture
    assert db.get_parent_tax_id(4) == 2, 'node with tax_id 4 has parent with tax_id 2'
    assert db.get_parent_tax_id(1) == None, 'the root node has no parent'
    
def test_get_tax_ids_path_root_to_node(setup_db_with_fixture):
    db = setup_db_with_fixture
    path = db.get_tax_ids_path_root_to_node(4)
    assert isinstance(path, list) 
    assert path == [1, 2, 4], 'correct path of tax_ids from root (1) to the requested node (4)'
    
def test_get_all_tax_ids_paths_root_to_leaf(setup_db_with_fixture):
    db = setup_db_with_fixture
    data = db.get_all_tax_ids_paths_root_to_leaf()
    assert isinstance(data, list), 'returns a list'
    assert len(data) >0, 'list is not empty'
    assert isinstance(data[0], list), 'elements are lists'
    assert [1, 2, 4] in data, 'one of the paths from root to leaf is 1->2->4'
    assert [1,3] in data , 'one of the paths from root to leaf is 1->3'
    assert len(data) ==2 , 'there are 2 paths from root to leaf nodes'
    
def test_get_flu_name_segment_data_dict(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    data = db.get_flu_name_segment_data_dict()
    assert isinstance(data, dict), 'returns a dict'
    assert len(data.keys())>0, 'we got some flu sequence data'
    expected_flu_names =set([
        'A/California/07/2009(H1N1)',
        'A/New York/392/2004(H3N2)',
        'A/Puerto Rico/8/1934(H1N1)',
        'B/Lee/1940'
    ])
    assert set(data.keys()) == expected_flu_names, 'the data contains the four expected flu names'
    assert set(data['A/Puerto Rico/8/1934(H1N1)'].keys())==set([1,2,3,4,5,6,7,8]), '...there are keys for each of the 8 segments for this flu isolate'
    assert data['A/Puerto Rico/8/1934(H1N1)'][3]== 2233, '...and has the correct length recorded for segment 3'
    
def test_mark_as_not_included(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    
    # before we run the process, no records in the fixture sequences are marked as include=0
    stmt='SELECT * FROM sequences WHERE include=0'
    rows = db._cur.execute(stmt).fetchall()
    assert len(rows) == 0 , 'there are no sequences marked as not included before we begin'
    ids=[1,5,8,10]
    db.mark_as_not_included(ids=ids)
    rows = db._cur.execute(stmt).fetchall()
    assert len(rows) == 4 , 'after marking as not included, we now retrieve 4 rows of sequences with included=0'

def test_retrieve_unnamed_unsegmented_flu(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    ids = db.retrieve_unnamed_unsegmented_flu()
    assert len(ids)==1, 'a single sequence in the fixture data is an "unnamed flu", ie no proper isolate name but identified as flu (kraken:taxid|518987|NC_002204.1 Influenza B virus RNA 1, complete sequence)'
    
def test_retrieve_sequence_ids_by_flu_name(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    ids = db.retrieve_sequence_ids_by_flu_name('A/Puerto Rico/8/1934(H1N1)')
    assert isinstance(ids, list), 'returns a list'
    assert len(ids) == 8 , 'there are 8 records with this name'
    
def test_max_tax_id(setup_db_with_fixture):
    db = setup_db_with_fixture
    assert db.max_tax_id() == 4,'the maximum tax_id in the fixtures nodes table is 4'
    
def test_retrieve_tax_id_by_node_scientific_name(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    assert db.retrieve_tax_id_by_node_scientific_name('Influenza A virus') == 11320,'correct tax_id found for fluA'
