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
        percent_n = 0,
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
    


def test_prune_db(setup_db):
    db = setup_db
    rows = db._cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='acc2taxids'").fetchall()
    assert len(rows)==1, 'the table acc2taxids exists before pruning'
    assert db.prune_db(), 'method returns True'
    rows = db._cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='acc2taxids'").fetchall()
    assert not rows, 'the table acc2taxids no longer exists after pruning'

def test_retrieve_flu_a_wo_subtype_ids(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    ids = db.retrieve_ids_flu_a_wo_subtype()
    assert not ids,'there are no cases of flu A without a subtype in the fixtures'
    
    # delete H subtype from one flu A genome in fixtures and re-retrieve
    stmt="""
    UPDATE sequences
    SET flu_a_h_subtype = NULL
    WHERE ncbi_acc = 'NC_002023.1'
    """
    db._cur.execute(stmt)
    db._con.commit()
    ids = db.retrieve_ids_flu_a_wo_subtype()
    assert len(ids) == 1,'having deleted the H subtype of one flu sequence, the method now retrieves one sequences.id'
    assert ids[0] == 1, 'the correct sequences.id is returned'
    
def test_sequences_category_exists(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    id=1
    label='test label'
    assert not db.sequences_category_exists(label), 'there are no matching sequences with this label in the fixtures'
    db._cur.execute("UPDATE sequences SET category= ? WHERE id= ?",[label, id])
    db._con.commit()
    assert db.sequences_category_exists(label), 'after the DB update, a sequence with the label exists in the DB and the method returns True'

def test_get_seq_ids_by_category_and_seq_lt(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    ids = db.get_seq_ids_by_category_and_seq_lt(category='test label', seq_len_lt=1000 )
    assert not ids, 'before we make changes to the fixtures, no sequences match the criteria'
    
    update_ids=[3,5,7]
    stmt="UPDATE sequences SET category = 'test label', seq_length = 900 WHERE id = ?"
    for id in update_ids:
        db._cur.execute(stmt,[id])
    db._con.commit()
    ids = db.get_seq_ids_by_category_and_seq_lt(category='test label', seq_len_lt=1000 )
    assert sorted(ids) == sorted(update_ids), 'after setting cateogry and seq length in three sequences to match filter, the three ids are returned'
    
def test_get_children_tax_ids(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    child_tax_ids = db.get_children_tax_ids(tax_id= 11250)
    assert sorted(child_tax_ids)==[208893,208895,410233], 'found the expected child tax_ids for taxonomy node 11250'
    
def test_get_sequence_ids_linked_to_taxon(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    seq_ids = db.get_sequence_ids_linked_to_taxon(tax_id= 2955291, include_children= True)
    assert not seq_ids, 'before we make changes to the fixtures, no sequences are linked to taxon 2955291 or any node in the sub-tree rooted at this taxon'
    
    # The relationships between the taxa in the fixtures are:                    
    #           2955291         
    #             │            
    #             │            
    #         ┌──11320─────┐    
    #         ▼            ▼    
    # ┌───114727─┐     119210  
    # │          │        │    
    # ▼          ▼        ▼    
    # 641809    211044    335341 
    #
    # We are linking sequences as follows:
    # tax_id    sequences.id
    # 11320     1
    # 114727    2
    # 211044    3
    # 641809    4
    # 335341    6
    # None of these are actually sequences of those taxa but this doesn't matter for the purpose of the test
    
    update_data = [
        [11320,1],
        [114727,2],
        [211044,3],
        [641809,4],
        [335341,6]
    ]
    db._cur.executemany("UPDATE sequences SET tax_id = ? WHERE id = ?", update_data)
    db._con.commit()
    
    # querying for sequence ids from taxon 2955291 "downwards" should now return the 5 sequence IDs from the update data
    seq_ids = db.get_sequence_ids_linked_to_taxon(tax_id= 2955291, include_children= True)
    assert sorted(seq_ids)==[1,2,3,4,6], 'having now linked 5 sequences to taxa in the sub-tree from 2955291, we return all 5 sequences.id'
    
    # repeat, but exclude taxon id 119210 and its children, which should result in sequence ID 6 not 
    # being in the resultset because it is linked to a child of 119210
    seq_ids = db.get_sequence_ids_linked_to_taxon(tax_id= 2955291, include_children= True, skip_tax_ids=[119210])
    assert sorted(seq_ids)==[1,2,3,4], 'having now linked 5 sequences to taxa in the sub-tree from 2955291, we return all 5 sequences.id'

def test_get_seq_ids_and_fasta_headers_by_category(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    data = db.get_seq_ids_and_fasta_headers_by_category(category='test label')
    assert not data, 'before we make changes to the fixtures, no sequences match the criteria'
    
    update_ids=[3,5,7]
    stmt="UPDATE sequences SET category = 'test label', fasta_header='test header' WHERE id = ?"
    for id in update_ids:
        db._cur.execute(stmt,[id])
    db._con.commit()
    data = db.get_seq_ids_and_fasta_headers_by_category(category='test label')
    assert isinstance(data, list), 'a list is returned'
    assert len(data) == 3, 'there are three items of data'
    assert {'id': 3, 'fasta_header': 'test header'} in data, 'returned data contains an expected sequence ID and FASTA header'
    
def test_get_duplicate_sequence_data(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    rows = db.get_duplicate_sequence_data()
    assert not rows, 'there are no duplicated sequences in the fixture set'
    
    # create two duplicates by inserting sequence data
    update_stmt = """
        INSERT INTO sequences (id,tax_id,fasta_header,dna_sequence,percent_n,seq_length,segment_number,ncbi_acc,flu_name,flu_type,flu_a_h_subtype,flu_a_n_subtype,include,is_flu,category,original_tax_id) VALUES
            (200,NULL,'kraken:taxid|211044|NC_002021.1 Influenza A virus (A/Puerto Rico/8/1934(H1N1)) segment 2, complete sequence','AGCGAAAGCAGGCAAACCATTTGAATGGATGTCAATCCGACCTTACTTTTCTTAAAAGTGCCAGCACAAAATGCTATAAGCACAACTTTCCCTTATACCGGAGACCCTCCTTACAGCCATGGGACAGGAACAGGATACACCATGGATACTGTCAACAGGACACATCAGTACTCAGAAAAGGCAAGATGGACAACAAACACCGAAACTGGAGCACCGCAACTCAACCCGATTGATGGGCCACTGCCAGAAGACAATGAACCAAGTGGTTATGCCCAAACAGATTGTGTATTGGAAGCAATGGCTTTCCTTGAGGAATCCCATCCTGGTATTTTTGAAAACTCGTGTATTGAAACGATGGAGGTTGTTCAGCAAACACGAGTAGACAAGCTGACACAAGGCCGACAGACCTATGACTGGACTTTAAATAGAAACCAGCCTGCTGCAACAGCATTGGCCAACACAATAGAAGTGTTCAGATCAAATGGCCTCACGGCCAATGAGTCTGGAAGGCTCATAGACTTCCTTAAGGATGTAATGGAGTCAATGAAAAAAGAAGAAATGGGGATCACAACTCATTTTCAGAGAAAGAGACGGGTGAGAGACAATATGACTAAGAAAATGATAACACAGAGAACAATAGGTAAAAGGAAACAGAGATTGAACAAAAGGAGTTATCTAATTAGAGCATTGACCCTGAACACAATGACCAAAGATGCTGAGAGAGGGAAGCTAAAACGGAGAGCAATTGCAACCCCAGGGATGCAAATAAGGGGGTTTGTATACTTTGTTGAGACACTGGCAAGGAGTATATGTGAGAAACTTGAACAATCAGGGTTGCCAGTTGGAGGCAATGAGAAGAAAGCAAAGTTGGCAAATGTTGTAAGGAAGATGATGACCAATTCTCAGGACACCGAACTTTCTTTGACCATCACTGGAGATAACACCAAATGGAACGAAAATCAGAATCCTCGGATGTTTTTGGCCATGATCACATATATGACCAGAAATCAGCCCGAATGGTTCAGAAATGTTCTAAGTATTGCTCCAATAATGTTCTCAAACAAAATGGCGAGACTGGGAAAAGGGTATATGTTTGAGAGCAAGAGTATGAAACTTAGAACTCAAATACCTGCAGAAATGCTAGCAAGCATTGATTTGAAATATTTCAATGATTCAACAAGAAAGAAGATTGAAAAAATCCGACCGCTCTTAATAGAGGGGACTGCATCATTGAGCCCTGGAATGATGATGGGCATGTTCAATATGTTAAGCACTGTATTAGGCGTCTCCATCCTGAATCTTGGACAAAAGAGATACACCAAGACTACTTACTGGTGGGATGGTCTTCAATCCTCTGACGATTTTGCTCTGATTGTGAATGCACCCAATCATGAAGGGATTCAAGCCGGAGTCGACAGGTTTTATCGAACCTGTAAGCTACATGGAATCAATATGAGCAAGAAAAAGTCTTACATAAACAGAACAGGTACATTTGAATTCACAAGTTTTTTCTATCGTTATGGGTTTGTTGCCAATTTCAGCATGGAGCTTCCCAGTTTTGGTGTGTCTGGGAGCAACGAGTCAGCGGACATGAGTATTGGAGTTACTGTCATCAAAAACAATATGATAAACAATGATCTTGGTCCAGCAACAGCTCAAATGGCCCTTCAGTTGTTCATCAAAGATTACAGGTACACGTACCGATGCCATAGAGGTGACACACAAATACAAACCCGAAGATCATTTGAAATAAAGAAACTGTGGGAGCAAACCCGTTCCAAAGCTGGACTGCTGGTCTCCGACGGAGGCCCAAATTTATACAACATTAGAAATCTCCACATTCCTGAAGTCTGCCTAAAATGGGAATTGATGGATGAGGATTACCAGGGGCGTTTATGCAACCCACTGAACCCATTTGTCAGCCATAAAGAAATTGAATCAATGAACAATGCAGTGATGATGCCAGCACATGGTCCAGCCAAAAACATGGAGTATGATGCTGTTGCAACAACACACTCCTGGATCCCCAAAAGAAATCGATCCATCTTGAATACAAGTCAAAGAGGAGTACTTGAAGATGAACAAATGTACCAAAGGTGCTGCAATTTATTTGAAAAATTCTTCCCCAGCAGTTCATACAGAAGACCAGTCGGGATATCCAGTATGGTGGAGGCTATGGTTTCCAGAGCCCGAATTGATGCACGGATTGATTTCGAATCTGGAAGGATAAAGAAAGAAGAGTTCACTGAGATCATGAAGATCTGTTCCACCATTGAAGAGCTCAGACGGCAAAAATAGTGAATTTAGCTTGTCCTTCATGAAAAAATGCCTTGTTCCTACT',0,2341,2,'NC_002021.1','A/Puerto Rico/8/1934(H1N1)','A',1,1,1,1,NULL,211044),
            (300,NULL,'kraken:taxid|211044|NC_002022.1 Influenza A virus (A/Puerto Rico/8/1934(H1N1)) segment 3, complete sequence','AGCGAAAGCAGGTACTGATCCAAAATGGAAGATTTTGTGCGACAATGCTTCAATCCGATGATTGTCGAGCTTGCGGAAAAAACAATGAAAGAGTATGGGGAGGACCTGAAAATCGAAACAAACAAATTTGCAGCAATATGCACTCACTTGGAAGTATGCTTCATGTATTCAGATTTCCACTTCATCAATGAGCAAGGCGAGTCAATAATCGTAGAACTTGGTGATCCTAATGCACTTTTGAAGCACAGATTTGAAATAATCGAGGGAAGAGATCGCACAATGGCCTGGACAGTAGTAAACAGTATTTGCAACACTACAGGGGCTGAGAAACCAAAGTTTCTACCAGATTTGTATGATTACAAGGAAAATAGATTCATCGAAATTGGAGTAACAAGGAGAGAAGTTCACATATACTATCTGGAAAAGGCCAATAAAATTAAATCTGAGAAAACACACATCCACATTTTCTCGTTCACTGGGGAAGAAATGGCCACAAAGGCCGACTACACTCTCGATGAAGAAAGCAGGGCTAGGATCAAAACCAGGCTATTCACCATAAGACAAGAAATGGCCAGCAGAGGCCTCTGGGATTCCTTTCGTCAGTCCGAGAGAGGAGAAGAGACAATTGAAGAAAGGTTTGAAATCACAGGAACAATGCGCAAGCTTGCCGACCAAAGTCTCCCGCCGAACTTCTCCAGCCTTGAAAATTTTAGAGCCTATGTGGATGGATTCGAACCGAACGGCTACATTGAGGGCAAGCTGTCTCAAATGTCCAAAGAAGTAAATGCTAGAATTGAACCTTTTTTGAAAACAACACCACGACCACTTAGACTTCCGAATGGGCCTCCCTGTTCTCAGCGGTCCAAATTCCTGCTGATGGATGCCTTAAAATTAAGCATTGAGGACCCAAGTCATGAAGGAGAGGGAATACCGCTATATGATGCAATCAAATGCATGAGAACATTCTTTGGATGGAAGGAACCCAATGTTGTTAAACCACACGAAAAGGGAATAAATCCAAATTATCTTCTGTCATGGAAGCAAGTACTGGCAGAACTGCAGGACATTGAGAATGAGGAGAAAATTCCAAAGACTAAAAATATGxxxxxxxCAAGTCAGCTAAAGTGGGCACTTGGTGAGAACATGGCACCAGAAAAGGTAGACTTTGACGACTGTAAAGATGTAGGTGATTTGAAGCAATATGATAGTGATGAACCAGAATTGAGGTCGCTTGCAAGTTGGATTCAGAATGAGTTCAACAAGGCATGCGAACTGACAGATTCAAGCTGGATAGAGCTTGATGAGATTGGAGAAGATGTGGCTCCAATTGAACACATTGCAAGCATGAGAAGGAATTATTTCACATCAGAGGTGTCTCACTGCAGAGCCACAGAATACATAATGAAGGGGGTGTACATCAATACTGCCTTACTTAATGCATCTTGTGCAGCAATGGATGATTTCCAATTAATTCCAATGATAAGCAAGTGTAGAACTAAGGAGGGAAGGCGAAAGACCAACTTGTATGGTTTCATCATAAAAGGAAGATCCCACTTAAGGAATGACACCGACGTGGTAAACTTTGTGAGCATGGAGTTTTCTCTCACTGACCCAAGACTTGAACCACACAAATGGGAGAAGTACTGTGTTCTTGAGATAGGAGATATGCTTCTAAGAAGTGCCATAGGCCAGGTTTCAAGGCCCATGTTCTTGTATGTGAGGACAAATGGAACCTCAAAAATTAAAATGAAATGGGGAATGGAGATGAGGCGTTGTCTCCTCCAGTCACTTCAACAAATTGAGAGTATGATTGAAGCTGAGTCCTCTGTCAAAGAGAAAGACATGACCAAAGAGTTCTTTGAGAACAAATCAGAAACATGGCCCATTGGAGAGTCTCCCAAAGGAGTGGAGGAAAGTTCCATTGGGAAGGTCTGCAGGACTTTATTAGCAAAGTCGGTATTTAACAGCTTGTATGCATCTCCACAACTAGAAGGATTTTCAGCTGAATCAAGAAAACTGCTTCTTATCGTTCAGGCTCTTAGGGACAATCTGGAACCTGGGACCTTTGATCTTGGGGGGCTATATGAAGCAATTGAGGAGTGCCTAATTAATGATCCCTGGGTTTTGCTTAATGCTTCTTGGTTCAACTCCTTCCTTACACATGCATTGAGTTAGTTGTGGCAGTGCTACTATTTGCTATCCATACTGTCCAAAAAAGTACCTTGTTTCTACT',0,2233,3,'NC_002022.1','A/Puerto Rico/8/1934(H1N1)','A',1,1,1,1,NULL,211044)
    """
    db._cur.executescript(update_stmt)
    db._con.commit()
    
    rows = db.get_duplicate_sequence_data()
    assert len(rows)==2, 'having inserted duplicates for 2 fixture sequences, we now retrieve 2 rows of duplicate results'
    assert rows[0]['min_id'] == 2,'first duplicate group has correct minimum id (the original record)'
    assert rows[0]['count'] == 2,'first duplicate group has correct count of duplicates in the group'
    assert rows[0]['ids'] == '2,200','first duplicate group has correct list of member ids'

def test_get_sequence_ids_percent_n_filter(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    seq_ids = db.get_sequence_ids_percent_n_filter(10)
    assert seq_ids == [] ,'there are no sequences with >10%N in the fixtures'
    
    # Set percent_n to values around the threshold for three new sequences.  
    # NOTE that the precent_n is not calculated from the sequence here, just provided as values to the INSERT.  
    # In a real DB, the percent_n value is calculated on ingest by the fasta loader
    insert_stmt = """
        INSERT INTO sequences (id,tax_id,fasta_header,dna_sequence,percent_n,seq_length,segment_number,ncbi_acc,flu_name,flu_type,flu_a_h_subtype,flu_a_n_subtype,include,is_flu,category,original_tax_id) VALUES
            (1000,NULL,'a sequence with 5% Ns','AAA',5,3,2,'NC_002021.1',NULL,NULL,NULL,NULL,1,0,NULL,NULL),
            (1001,NULL,'a sequence with 10% Ns','AAA',10,3,2,'NC_002021.1',NULL,NULL,NULL,NULL,1,0,NULL,NULL),
            (1002,NULL,'a sequence with 15% Ns','AAA',15,3,2,'NC_002021.1',NULL,NULL,NULL,NULL,1,0,NULL,NULL)
    """
    db._cur.executescript(insert_stmt)
    db._con.commit()
    
    seq_ids = db.get_sequence_ids_percent_n_filter(10)
    assert seq_ids == [1002] ,'there is now one sequence in the DB with >10% N and the ID is correct'
    
def test_tax_id_exists(setup_db_with_fixture):
    db = setup_db_with_fixture
    assert db.tax_id_exists(1), 'tax_id 1 exists in the DB'
    assert not db.tax_id_exists(1000000), 'tax_id 1000000 does not exist in the DB'
