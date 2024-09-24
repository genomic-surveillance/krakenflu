import pytest
import os.path
from importlib_resources import files

from kraken_flu.src.db import Db, BulkInsertBuffer

# set this to True to make the tests print all executed SQL or False to stop that
PRINT_SQL_TRACE=True

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

def test_add_taxon(setup_db_with_fixture):
    db = setup_db_with_fixture
    stmt="""
    SELECT 
        taxonomy_nodes.parent_tax_id AS parent_tax_id,
        taxonomy_names.name AS name,
        taxonomy_names.name_class AS name_class
    FROM taxonomy_nodes
    INNER JOIN taxonomy_names ON(taxonomy_nodes.tax_id = taxonomy_names.tax_id)
    WHERE taxonomy_nodes.tax_id = ?
    """
    tax_id = 1000
    parent_tax_id = 1
    name='a new taxon'
    rows = db._cur.execute(stmt,[tax_id]).fetchall()
    assert not rows, 'before we insert the new taxon, it does not exist'
    db.add_taxon(tax_id= tax_id, parent_tax_id= parent_tax_id, name= name) 
    rows = db._cur.execute(stmt,[tax_id]).fetchall()
    assert len(rows)==1, 'now the new taxon exists'
    assert rows[0]['name'] == name, '...it has the correct name'
    assert rows[0]['parent_tax_id'] == parent_tax_id, '...it has the correct parent'
    assert rows[0]['name_class'] == 'scientific name', '...it has the hardcoded name_class scientific name'
    
def test_retrieve_all_flu_sequences(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    rows = db.retrieve_all_flu_sequences()
    assert isinstance(rows, list) 
    assert len(rows)==32, 'there are 32 records of sequences with the is_flu flag set'
    
def test_set_tax_id_for_sequence(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    stmt="SELECT tax_id FROM sequences WHERE id=1"
    row = db._cur.execute(stmt).fetchone()
    assert row
    assert not row['tax_id'], 'before we start, sequences.id=1 has no tax_id set'
    db.set_tax_id_for_sequence(id=1, tax_id=999)
    row = db._cur.execute(stmt).fetchone()
    assert row['tax_id']==999, 'after the update, the record has the correct tax_id set'
    
def test_all_sequences_iterator(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    it = db.all_sequences_iterator()
    row=next(it)
    assert row, 'can retrieve a row at a time from iterator'
    assert row['fasta_header'], 'the row has a fasta-Header attribute'
    
    # get the remaining rows
    rows=[]
    for row in it:
        rows.append(row)
    assert len(rows) == 34, 'there are 35 sequences in the fixtures of which 34 remain after the first row has already been retrieved before'
    
    # start again and fetch all 35 rows into the list
    it = db.all_sequences_iterator()

    rows=[]
    for row in it:
        rows.append(row)
    assert len(rows) == 35, 'retrieved all 35 sequences in the fixtures'
    
def test_get_field_names(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    expected_fields = [
        'id',
        'tax_id', 
        'name',
        'name_class', 
        'unique_name'
    ]
    assert db.get_field_names('taxonomy_names') == expected_fields, 'retrieving the expected field names for a table'

def test_bulk_insert(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    
    rows = db._cur.execute("SELECT * FROM taxonomy_names WHERE name IN(?,?)",['bulk1','bulk2']).fetchall()
    assert not rows, 'before we do the bulk insert, the rows do not exist in the table'
    
    table_name= 'taxonomy_names'
    field_names= [
        'tax_id', 
        'name',
        'name_class', 
        'unique_name'
    ]
    data= [
        [ 100001, 'bulk1', 'scientific name', None],
        [ 100002, 'bulk2', 'scientific name', None]
    ]
    db.bulk_insert(table_name= table_name, field_names= field_names, field_data= data)
    
    rows = db._cur.execute("SELECT * FROM taxonomy_names WHERE name IN(?,?)",['bulk1','bulk2']).fetchall()
    assert len(rows) == 2, 'after we do the bulk insert, the 2 new rows exist in the table'
    
def test_bulk_insert_buffer(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    
    rows = db._cur.execute("SELECT * FROM taxonomy_names WHERE name IN(?,?)",['bulk1','bulk2']).fetchall()
    assert not rows, 'before we do the bulk insert, the rows do not exist in the table'

    data= [
        { 'tax_id': 100001, 'name': 'bulk1', 'name_class': 'scientific name', 'unique_name': None},
        { 'tax_id': 100002, 'name': 'bulk2', 'name_class': 'scientific name', 'unique_name': None}
    ]
    
    with db.bulk_insert_buffer('taxonomy_names') as b:
        for row in data:
            n_inserted = b.add_row(row)
            assert n_inserted==0, 'buffer size is 5000, so every add_row call just buffers (will flush to DB when losing context)'
    
    rows = db._cur.execute("SELECT * FROM taxonomy_names WHERE name IN(?,?)",['bulk1','bulk2']).fetchall()
    assert len(rows) == 2, 'after we do the bulk insert, the 2 new rows exist in the table'
    
def test_bulk_insert_buffer_force_single(setup_db_with_real_world_fixture):
    # same as test_bulk_insert_buffer but set the buffer size to single row, so that each row will trigger one INSERT, 
    # thus testing the ability to buffer batches of data and flush to DB
    
    db = setup_db_with_real_world_fixture
    rows = db._cur.execute("SELECT * FROM taxonomy_names WHERE name IN(?,?)",['bulk1','bulk2']).fetchall()
    assert not rows, 'before we do the bulk insert, the rows do not exist in the table'

    data= [
        { 'tax_id': 100001, 'name': 'bulk1', 'name_class': 'scientific name', 'unique_name': None},
        { 'tax_id': 100002, 'name': 'bulk2', 'name_class': 'scientific name', 'unique_name': None}
    ]
    
    with db.bulk_insert_buffer(table_name='taxonomy_names', buffer_size=1) as b:
        for row in data:
            n_inserted = b.add_row(row)
            assert n_inserted == 1, 'every row added triggers a flush to DB'
    
    rows = db._cur.execute("SELECT * FROM taxonomy_names WHERE name IN(?,?)",['bulk1','bulk2']).fetchall()
    assert len(rows) == 2, 'after we do the bulk insert, the 2 new rows exist in the table'
    
def test_bulk_update(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    
    # get some sequences records from the fixtures
    id_field_values = [4,5,7]
    rows = db._cur.execute("SELECT id, tax_id, mod_fasta_header FROM sequences WHERE id IN(?,?,?) ORDER BY id",id_field_values).fetchall()
    assert len(rows)==3, '3 rows of sequecnces data have been retrieved'
    assert [x['tax_id'] for x in rows ] == [None, None, None], 'before the bulk update, none of the records has a tax_id'
    assert [x['mod_fasta_header'] for x in rows ] == [None, None, None], 'before the bulk update, none of the records has a mod_fasta_header'

    field_data = [ [200,'new header 1'], [300,'new header 2'],[456,'another new header']]

    # do the bulk update
    db.bulk_update(
        table_name= 'sequences',
        update_fields= ['tax_id', 'mod_fasta_header'],
        id_field = 'id',
        field_data= field_data,
        id_field_values= id_field_values
    )
    
    # query for the records again, they should now have been updated
    rows = db._cur.execute("SELECT id, tax_id, mod_fasta_header FROM sequences WHERE id IN(?,?,?) ORDER BY id",id_field_values).fetchall()
    assert len(rows)==3, '3 rows of sequecnces data have been retrieved'
    assert [x['tax_id'] for x in rows ] == [200, 300, 456], 'after the bulk update, the records have the expected new tax_id values'
    assert [x['mod_fasta_header'] for x in rows ] == ['new header 1', 'new header 2', 'another new header'], 'before the bulk update, none of the records has a mod_fasta_header'


def test_bulk_update_buffer(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    id_field_values = [4,5,7]
    rows = db._cur.execute("SELECT id, tax_id, mod_fasta_header FROM sequences WHERE id IN(?,?,?) ORDER BY id",id_field_values).fetchall()
    assert len(rows)==3, '3 rows of sequecnces data have been retrieved'
    assert [x['tax_id'] for x in rows ] == [None, None, None], 'before the bulk update, none of the records has a tax_id'
    assert [x['mod_fasta_header'] for x in rows ] == [None, None, None], 'before the bulk update, none of the records has a mod_fasta_header'

    field_data= [
        { 'id': 4, 'tax_id': 200, 'mod_fasta_header': 'new header 1'},
        { 'id': 5, 'tax_id': 300, 'mod_fasta_header': 'new header 2'},
        { 'id': 7, 'tax_id': 456, 'mod_fasta_header': 'another new header'}
    ]
    
    # run the buffered updates
    # the buffer size is 5000, so the three updates will all be carried out when the context manager
    # loses focus while each add_row operation only adds to the data in memory and does not trigger a DB action
    with db.bulk_update_buffer(table_name='sequences', id_field='id', update_fields=['tax_id', 'mod_fasta_header'], buffer_size= 5000) as b:
        for row in field_data:
            n_updated = b.add_row(row)
            assert n_updated==0, 'buffer size is 5000, so every add_row call just buffers (will flush to DB when losing context)'
    
    # query for the records again, they should now have been updated
    rows = db._cur.execute("SELECT id, tax_id, mod_fasta_header FROM sequences WHERE id IN(?,?,?) ORDER BY id",id_field_values).fetchall()
    assert len(rows)==3, '3 rows of sequecnces data have been retrieved'
    assert [x['tax_id'] for x in rows ] == [200, 300, 456], 'after the bulk update, the records have the expected new tax_id values'
    assert [x['mod_fasta_header'] for x in rows ] == ['new header 1', 'new header 2', 'another new header'], 'before the bulk update, none of the records has a mod_fasta_header'

def test_bulk_update_buffer_force_single(setup_db_with_real_world_fixture):
    # same as test_bulk_update_buffer but set the buffer size to single row, so that each row will trigger one INSERT, 
    # thus testing the ability to buffer batches of data and flush to DB
    db = setup_db_with_real_world_fixture
    id_field_values = [4,5,7]
    rows = db._cur.execute("SELECT id, tax_id, mod_fasta_header FROM sequences WHERE id IN(?,?,?) ORDER BY id",id_field_values).fetchall()
    assert len(rows)==3, '3 rows of sequecnces data have been retrieved'
    assert [x['tax_id'] for x in rows ] == [None, None, None], 'before the bulk update, none of the records has a tax_id'
    assert [x['mod_fasta_header'] for x in rows ] == [None, None, None], 'before the bulk update, none of the records has a mod_fasta_header'

    field_data= [
        { 'id': 4, 'tax_id': 200, 'mod_fasta_header': 'new header 1'},
        { 'id': 5, 'tax_id': 300, 'mod_fasta_header': 'new header 2'},
        { 'id': 7, 'tax_id': 456, 'mod_fasta_header': 'another new header'}
    ]
    
    # run the buffered updates
    # the buffer size is 1, so this time, each of the add_row commands will trigger a DB update transaction
    with db.bulk_update_buffer(table_name='sequences', id_field='id', update_fields=['tax_id', 'mod_fasta_header'], buffer_size= 1) as b:
        for row in field_data:
            n_updated = b.add_row(row)
            assert n_updated==1, 'buffer size is 1, so every add_row call flushes 1 update to the DB'
    
    # query for the records again, they should now have been updated
    rows = db._cur.execute("SELECT id, tax_id, mod_fasta_header FROM sequences WHERE id IN(?,?,?) ORDER BY id",id_field_values).fetchall()
    assert len(rows)==3, '3 rows of sequecnces data have been retrieved'
    assert [x['tax_id'] for x in rows ] == [200, 300, 456], 'after the bulk update, the records have the expected new tax_id values'
    assert [x['mod_fasta_header'] for x in rows ] == ['new header 1', 'new header 2', 'another new header'], 'before the bulk update, none of the records has a mod_fasta_header'

def test_all_seq2taxid_iterator(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture

    acc2taxids_rows = db._cur.execute("SELECT accession,tax_id FROM acc2taxids").fetchall()
    assert len(acc2taxids_rows) == 4, '4 rows of acc2taxids in fixtures'

    it = db.all_seq2taxid_iterator()
    rows=[]
    for row in it:
        rows.append(row)
    # fixtures contains 4 rows of acc2taxids with two unique tax_ids, so in total we should link 
    # 2 sequences via accession ID to tax_id
    assert len(rows) == 2, 'all_seq2taxid_iterator retrieved all (2) rows linking sequences to tax IDs in acc2taxids table via NCBI accession'

    # check row data for the expected associations of sequences (by sequences.id) to tax_id
    # the first is for 'NC_001803.1 Respiratory syncytial virus, complete genome' which should be linked to tax_id 12814
    assert [x for x in rows if x['id']==9 and x['tax_id']==12814], 'the resultset contains an expected association of sequences.id 9 with tax_id 12814'

    # the second one is 'NC_002205.1 Influenza B virus (B/Lee/1940) segment 2', linked to tax_id 518987, sequences.id=21
    assert [x for x in rows if x['id']==21 and x['tax_id']==518987], 'the resultset contains an expected association of sequences.id 21 with tax_id 518987'
    
    # setting a tax_id for one of the two sequences linked above will make it drop from the 
    # resultset, which retrieves only unlinked (no tax_id yet) records by default
    db._cur.execute("UPDATE sequences SET tax_id=7878787 WHERE id=21")
    it = db.all_seq2taxid_iterator()
    rows=[]
    for row in it:
        rows.append(row)
    assert len(rows) == 1, 'after setting the tax_id field to a value for one of the previously returned sequences, it is no longer in the resultset'
    assert [x for x in rows if x['id']==9 and x['tax_id']==12814], 'the resultset still contains a linkage for sequences.id 9'
    assert not [x for x in rows if x['id']==21 and x['tax_id']==518987], 'sequences.id 21 no longer in the resultset because it already has a tax_id now'
    
def test_sequences_category_exists(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    id=1
    label='test label'
    assert not db.sequences_category_exists(label), 'there are no matching sequences with this label in the fixtures'
    db._cur.execute("UPDATE sequences SET category= ? WHERE id= ?",[label, id])
    db._con.commit()
    assert db.sequences_category_exists(label), 'after the DB update, a sequence with the label exists in the DB and the method returns True'
