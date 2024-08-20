import pytest
import os.path
from importlib_resources import files

from kraken_flu.src.kraken_db_builder import KrakenDbBuilder
from kraken_flu.src.db import Db

FIXTURE_DIR = files('kraken_flu.tests.fixtures')
TAX_DIR = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','taxonomy'))

SMALL_VIRUS_FILE = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','library','viral','library.fna'))
OTHER_SMALL_VIRUS_FILE = FIXTURE_DIR.joinpath(os.path.join('ncbi_data_not_kraken','library','viral','library.fna'))

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

def test_init_no_patht():
    kdb = KrakenDbBuilder()
    assert kdb
    assert isinstance(kdb._db, Db)
    assert kdb.db_path,'a temp file path has been assigned to build the DB'
    assert not kdb.db_ready(), 'the initial state of the DB is not ready'
    
def test_init_path( tmp_path ):
    db_path = tmp_path / 'kraken_flu.db'
    kdb = KrakenDbBuilder(db_path= db_path)
    assert kdb.db_path == db_path, 'providing a path results in the DB path being set to that path'    
    
def test_load_fasta():
    kdb = KrakenDbBuilder()
    assert kdb.fasta_files_loaded == [], 'before we upload a FASTA file, the list of loaded FASTA files is empty'
    kdb.load_fasta_file(file_path=SMALL_VIRUS_FILE)
    assert kdb.fasta_files_loaded == [SMALL_VIRUS_FILE ], 'after loading one file, the list of loaded FASTA contains the path of the single file loaded'

    # add a second FASTA file, set the category this time
    kdb.load_fasta_file(file_path=OTHER_SMALL_VIRUS_FILE, category='other')
    assert kdb.fasta_files_loaded == [SMALL_VIRUS_FILE, OTHER_SMALL_VIRUS_FILE], 'now there are two FASTA file path recorded as loaded into DB'

def test_load_taxonomy():
    kdb = KrakenDbBuilder()
    assert not kdb.taxonomy_loaded, 'taxonomy_loaded returns False before we load anything'
    kdb.load_taxonomy_files(taxonomy_dir=TAX_DIR)
    assert kdb.taxonomy_loaded, 'after loading, taxonomy_loaded returns True'
    
def test_db_ready():
    kdb = KrakenDbBuilder()
    assert not kdb.db_ready(), 'db_ready is False before we load the files into DB'
    kdb.load_taxonomy_files(taxonomy_dir=TAX_DIR)
    assert not kdb.db_ready(), 'db_ready is still False after loading just the taxonomy files into DB'
    kdb.load_fasta_file(file_path=SMALL_VIRUS_FILE)
    assert kdb.db_ready(), 'after also loading at least one FASTA file, db_ready is now True'
    
def test_filter_incomplete_flu(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    kdb = KrakenDbBuilder(db=db)
    not_included_stmt='SELECT * FROM sequences WHERE include=0'
    rows = db._cur.execute(not_included_stmt).fetchall()
    assert len(rows) == 0 , 'there are no sequences marked as not included before we begin'
    
    assert kdb.filter_incomplete_flu()
    
    # The data for flu B isolate B/Lee/1940 has a segment 1 with a name that cannot be parsed
    # For that reason, without any attempt to rescue this isolate, the method should now identify
    # 7 segments that do have a flu_name but this isolate is incomplete due to the missing segment 1
    rows = db._cur.execute(not_included_stmt).fetchall()
    assert len(rows) == 7 , 'due to one flu B seg1 not having a parsable flu isolate name, the 7 remaining segments are incomplete and hence marked not to be included'
    
    # Make another isolate incomplete by setting the segment length to something very small
    select_a_seg1_stmt = """
    SELECT id 
    FROM sequences
    WHERE flu_name = ?
    AND segment_number = 1
    """
    seq1_row = db._cur.execute(select_a_seg1_stmt,['A/Puerto Rico/8/1934(H1N1)']).fetchone()
    assert seq1_row
    update_stmt = """
    UPDATE sequences
    SET seq_length = 100
    WHERE id = ?
    """
    db._cur.execute(update_stmt,[seq1_row['id']])
    db._con.commit()

    # rerun the filter, it should now filter out all segments of this isolate because it no longer
    # has 8 full-length segments
    assert kdb.filter_incomplete_flu()

    rows = db._cur.execute(not_included_stmt).fetchall()
    assert len(rows) == 15 , 'setting seg1 one one isolate to a very small length results in all 8 segments of this isolate being excluded, bringing total exclude count to 15'
    
    # make another isolate incomplete but this time exclude it from the filter, so it should still pass
    isolate_name='A/California/07/2009(H1N1)'
    seq1_row = db._cur.execute(select_a_seg1_stmt,[isolate_name]).fetchone()
    assert seq1_row
    db._cur.execute(update_stmt,[seq1_row['id']])
    db._con.commit()
    
    assert kdb.filter_incomplete_flu(filter_except_patterns=[isolate_name])
    rows = db._cur.execute(not_included_stmt).fetchall()
    assert len(rows) == 15 , 'making another isolate incomplete but also including it in the exempt names should leave the isolate unfiltered and the exclude count unchanged'
    
def test_next_new_tax_id(setup_db_with_fixture):
    db = setup_db_with_fixture
    kdb = KrakenDbBuilder(db=db)
    assert kdb.next_new_tax_id() == db.max_tax_id() + 1, 'the initial next new tax_id is the maximum tax_id from the DB +1'
    assert kdb.next_new_tax_id() == db.max_tax_id() + 2, 'requesting a new tax id again returns an incremented id'

def test_create_segmented_flu_taxonomy_nodes(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    kdb = KrakenDbBuilder(db=db)
    stmt = """
    SELECT 
        taxonomy_nodes.tax_id AS tax_id,
        taxonomy_nodes.parent_tax_id AS parent_tax_id,
        taxonomy_names.name AS name,
        taxonomy_names.name_class AS name_class
    FROM taxonomy_nodes
    INNER JOIN taxonomy_names ON(taxonomy_nodes.tax_id = taxonomy_names.tax_id)
    WHERE taxonomy_names.name = ? and taxonomy_names.name_class = 'scientific name'
    """
    rows = db._cur.execute(stmt,['Influenza A virus']).fetchall()
    assert len(rows) == 1, 'a node exists for "Influenza Virus A" before we start'
    rows = db._cur.execute(stmt,['Influenza A segment 1']).fetchall()
    assert not rows, 'no node exists for "Influenza A segment 1" before we start'
    
    new_node_ids = kdb.create_segmented_flu_taxonomy_nodes()
    assert isinstance(new_node_ids, dict), 'create_segmented_flu_taxonomy_node returns a dict'
    assert isinstance(new_node_ids['A'][None][1],int), 'the data structure returned assigns an int to a combination of virus type subtype and segment number'
    type_A_seg1_tax_id = new_node_ids['A'][None][1]
    assert type_A_seg1_tax_id >0 , 'the type A segment 1 new tax_id is a non-zero integer'
    
    rows = db._cur.execute(stmt,['Influenza A segment 1']).fetchall()
    assert len(rows), 'a node exists now for "Influenza A segment 1"'
    assert rows[0]['tax_id'] == new_node_ids['A'][None][1], 'the new node tax_id has been correctly recorded in the returned datastrcuture'
    
    rows = db._cur.execute(stmt,['Influenza A segment 4']).fetchall()
    assert len(rows), 'a node exists now for "Influenza A segment 4"'
    infA_seg_4_tax_id = new_node_ids['A'][None][4]
    assert infA_seg_4_tax_id >0 , 'got a tax_id for the new node "Influenza A segment 4"'
    
    rows = db._cur.execute(stmt,['Influenza A H2 segment 4']).fetchall()
    assert len(rows), 'a node exists now for "Influenza H2 subtype segment 4"'
    infA_H2_seg_4_tax_id = new_node_ids['A']['H2'][4]
    assert infA_H2_seg_4_tax_id >0 , 'got a tax_id for the new node "Influenza A H2 segment 4"'
    infA_H2_seg_4_parent_tax_id=rows[0]['parent_tax_id']
    assert infA_H2_seg_4_parent_tax_id>0, 'got a parent_tax_id for  "Influenza A H2 segment 4"'

    assert infA_H2_seg_4_parent_tax_id == infA_seg_4_tax_id ,'the parent of "Influenza A H2 segment 4" is  "Influenza A segment 4"'

    rows = db._cur.execute(stmt,['Influenza A N11 segment 6']).fetchall()
    assert len(rows), 'a node exists now for "Influenza N11 subtype segment 6"'
    
    rows = db._cur.execute(stmt,['Influenza B segment 4']).fetchall()
    assert len(rows), 'a node exists now for "Influenza B segment 4"'
    
    rows = db._cur.execute(stmt,['Influenza B H2 segment 4']).fetchall()
    assert not rows, 'for influenza B, no H2 subtype nodes under segment 4 have been created"'
    
    rows = db._cur.execute(stmt,['Influenza C segment 1']).fetchall()
    assert not rows, 'no Influenza C parent node exists in the fixtures, hence no Influenza C new taxonomy a node exists now for "Influenza H2 subtype segment 4"'
    
def test_assign_flu_taxonomy_nodes(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    kdb = KrakenDbBuilder(db=db)
    stmt1 = """
        SELECT 
            sequences.tax_id AS seq_tax_id,
            sequences.mod_fasta_header AS mod_fasta_header,
            taxonomy_nodes.tax_id AS node_tax_id,
            taxonomy_nodes.parent_tax_id AS seq_parent_tax_id
        FROM sequences
        LEFT OUTER JOIN taxonomy_nodes ON (sequences.tax_id = taxonomy_nodes.tax_id)
        WHERE flu_name = ? AND segment_number = ?
    """
    stmt2 = """
        SELECT 
            taxonomy_nodes.tax_id AS tax_id,
            taxonomy_names.name
        FROM taxonomy_nodes
        INNER JOIN taxonomy_names ON(taxonomy_nodes.tax_id = taxonomy_names.tax_id)
        WHERE taxonomy_names.name = ?
    """
    
    fluA_puerto_rico_8_seq_rows = db._cur.execute(stmt1,["A/Puerto Rico/8/1934(H1N1)", 8]).fetchall()
    assert len(fluA_puerto_rico_8_seq_rows) == 1, 'found the sequences record for A/Puerto Rico/8/1934(H1N1) segment 8'
    assert not fluA_puerto_rico_8_seq_rows[0]['seq_tax_id'], '...before we start, no tax_id is assigned to this sequence'
    assert not fluA_puerto_rico_8_seq_rows[0]['node_tax_id'], '...no taxonomy_nodes record exists yet for this sequence'
    
    # run the tax_id assignment and re-query for the flu A sequence, which should now have an associated 
    # taxonomy node created, which should be a child node of one of the new flu A segment taxonomy nodes
    kdb.assign_flu_taxonomy_nodes()
    fluA_puerto_rico_8_rows = db._cur.execute(stmt1,["A/Puerto Rico/8/1934(H1N1)", 8]).fetchall()
    assert fluA_puerto_rico_8_rows[0]['seq_tax_id'], '...after running assign_flu_taxonomy_nodes, the record now has a tax_id assigned'
    assert fluA_puerto_rico_8_rows[0]['node_tax_id'], '...there is also a taxonomy_nodes record associated with this sequence now'
    assert fluA_puerto_rico_8_rows[0]['seq_parent_tax_id'], '...the associated taxonomy node has a parent_tax_id'
    assert fluA_puerto_rico_8_rows[0]['mod_fasta_header'] == 'A/Puerto Rico/8/1934(H1N1) segment 8', '...the sequences record has an alternative FASTA header assigned correctly'

    
    # retrieve the taxonomy node Inlfuenza A segment 8 and check that it is assigned as a parent 
    fluA_seg8_node = db._cur.execute(stmt2,['Influenza A segment 8']).fetchone()
    assert fluA_seg8_node, 'found flu A segment 8 node'
    assert fluA_puerto_rico_8_rows[0]['seq_parent_tax_id'] == fluA_seg8_node['tax_id'], '...the parent of this flu A segment 8 sequence is the new node "Influenza A segment 8"'
    
    # same for flu A segment 4, which should be the child of a "Influenza A H1 segment 4" node
    fluA_puerto_rico_4_rows = db._cur.execute(stmt1,["A/Puerto Rico/8/1934(H1N1)", 4]).fetchall()
    assert fluA_puerto_rico_4_rows[0]['seq_tax_id'], '...after running assign_flu_taxonomy_nodes, the record now has a tax_id assigned'
    assert fluA_puerto_rico_4_rows[0]['node_tax_id'], '...there is also a taxonomy_nodes record associated with this sequence now'
    assert fluA_puerto_rico_4_rows[0]['seq_parent_tax_id'], '...the associated taxonomy node has a parent_tax_id'
    
    # retrieve the taxonomy node Inlfuenza A segment 8 and check that it is assigned as a parent 
    fluA_h1_seg4_node = db._cur.execute(stmt2,['Influenza A H1 segment 4']).fetchone()
    assert fluA_h1_seg4_node, 'found flu A H1 segment 4 node'
    assert fluA_puerto_rico_4_rows[0]['seq_parent_tax_id'] == fluA_h1_seg4_node['tax_id'], '...the parent of this flu A segment 8 sequence is the new node "Influenza A segment 8"'
    
    
    # check a flu B segment, where there is no Hx segment 4, just Flu B segment 4 as the parent node
    fluB_4_rows = db._cur.execute(stmt1,["B/Lee/1940", 4]).fetchall()
    assert fluB_4_rows[0]['seq_tax_id'], 'Flu B/Lee/1940 has a tax_id assigned'
    assert fluB_4_rows[0]['node_tax_id'], '...there is also a taxonomy_nodes record associated with this sequence now'
    assert fluB_4_rows[0]['seq_parent_tax_id'], '...the associated taxonomy node has a parent_tax_id'
    fluB_seg4_node = db._cur.execute(stmt2,['Influenza B segment 4']).fetchone()
    assert fluA_seg8_node, 'found flu B segment 4 node'
    assert fluB_4_rows[0]['seq_parent_tax_id'] == fluB_seg4_node['tax_id'], '...the parent of this flu B segment 4 sequence is the new node "Influenza B segment 4"'
    