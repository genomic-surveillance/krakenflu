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
    
@pytest.fixture(scope='function')
def setup_db_with_rsv_fixture( setup_db ):
    db = setup_db
    fixture_dir = files('kraken_flu.tests.fixtures')
    db_fixture_file = fixture_dir.joinpath(os.path.join('rsv','rsv_db_fixture1.sql'))
    with open(db_fixture_file, 'r') as file:
        file_content = file.read()
    db._cur.executescript(file_content)
    yield db

# TODO need to clean up this duplication of fixture definitions
@pytest.fixture(scope='function')
def setup_db_with_rsv_fixture_for_filter( setup_db ):
    db = setup_db
    fixture_dir = files('kraken_flu.tests.fixtures')
    db_fixture_file = fixture_dir.joinpath(os.path.join('rsv','rsv_db_fixture2.sql'))
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
    kdb.load_taxonomy_files(taxonomy_dir=TAX_DIR, no_acc2taxid= True)
    assert kdb.taxonomy_loaded, 'after loading, taxonomy_loaded returns True'
    
def test_db_ready():
    kdb = KrakenDbBuilder()
    assert not kdb.db_ready(), 'db_ready is False before we load the files into DB'
    kdb.load_taxonomy_files(taxonomy_dir=TAX_DIR, no_acc2taxid= True)
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
    
def test_create_db_ready_dir(setup_db_with_real_world_fixture, tmp_path):
    db = setup_db_with_real_world_fixture
    kdb = KrakenDbBuilder(db=db)
    
    out_dir = tmp_path / "kdb_out_dir"
    exp_names_out_file = out_dir / 'taxonomy' / 'names.dmp'
    assert not exp_names_out_file.is_file(), 'before we start, the expected file names.dmp does not exist'
    
    exp_lib_out_file = out_dir / 'library' / 'library.fna'
    assert not exp_lib_out_file.is_file(), 'before we start, the expected file library.fna does not exist'
    
    kdb.create_db_ready_dir(path= out_dir)

    assert exp_names_out_file.is_file(), 'after running create_db_ready_dir, the expected file taxonomy/names.dmp exists'
    assert exp_lib_out_file.is_file(), 'after running create_db_ready_dir, the expected file library/library.fna exists'

def test_link_all_unlinked_sequences_to_taxonomy_nodes(setup_db_with_real_world_fixture, tmp_path):
    db = setup_db_with_real_world_fixture
    kdb = KrakenDbBuilder(db=db)
    
    stmt_base= """
        SELECT 
            sequences.id,
            sequences.fasta_header,
            sequences.tax_id,
            taxonomy_names.name AS taxonomy_name
        FROM sequences
        LEFT OUTER JOIN taxonomy_names ON(taxonomy_names.tax_id = sequences.tax_id AND taxonomy_names.name_class='scientific name')
    """
    stmt_tax_id_not_set = stmt_base + ' WHERE sequences.tax_id IS NULL'
    stmt_tax_id_set = stmt_base + ' WHERE sequences.tax_id IS NOT NULL'

    rows = db._cur.execute(stmt_tax_id_not_set).fetchall()
    assert len(rows) == 35, 'before running the linkage, all 35 sequences records have no tax_id set'
    
    # run the linkage on NCBI accession ID
    # only two sequences will be linked because the fixture data only contains acc2taxid records for two sequences:
    #   - 'NC_002205.1 Influenza B virus (B/Lee/1940) segment 2'
    #   - 'NC_001803.1 Respiratory syncytial virus, complete genome'
    kdb.link_all_unlinked_sequences_to_taxonomy_nodes()
    rows = db._cur.execute(stmt_tax_id_not_set).fetchall()
    assert len(rows) == 33, 'after running the linkage, 2 sequences have their tax_id set due to acc2taxid records in the fixtures while 33 sequences remain unlinked'
    
    rows = db._cur.execute(stmt_tax_id_set).fetchall()
    assert len(rows) == 2, 'check with reverse logic that 2 sequences now have a tax_id'
    assert [x for x in rows if x['id']==9 and x['tax_id']==12814], 'the resultset contains an expected association of sequences.id 9 with tax_id 12814'

def test_filter_flu_a_wo_subtype(setup_db_with_real_world_fixture, tmp_path):
    db = setup_db_with_real_world_fixture
    kdb = KrakenDbBuilder(db=db)
    
    # retrieve data for a flu A sequence which is not filtered out at this point
    rows = db._cur.execute("SELECT id, include FROM sequences WHERE ncbi_acc = 'NC_002023.1'").fetchall()
    assert len(rows) == 1, 'a flu A sequence can be retrieved from the DB by NCBI acc'
    assert rows[0]['include'] == 1, '... the sequence is set to be included'
    
    # run the filter - it should not affect any sequences at this point
    kdb.filter_flu_a_wo_subtype()
    rows = db._cur.execute("SELECT id, include FROM sequences WHERE ncbi_acc = 'NC_002023.1'").fetchall()
    assert rows[0]['include'] == 1, 'after running filter_flu_a_wo_subtype, this sequence is unaffected because it has H/N subtypes'

    # set the flu A H subtype of this sequence to NULL and re-run the filter
    db._cur.execute("UPDATE sequences SET flu_a_h_subtype = NULL WHERE ncbi_acc = 'NC_002023.1'")
    db._con.commit()
    kdb.filter_flu_a_wo_subtype()
    rows = db._cur.execute("SELECT id, include FROM sequences WHERE ncbi_acc = 'NC_002023.1'").fetchall()
    assert rows[0]['include'] == 0, 'having set H subtype to NULL, this sequence is now filtered out by filter_flu_a_wo_subtype'
    
def test__apply_size_filter_to_labelled_sequences(setup_db_with_real_world_fixture, tmp_path):
    db = setup_db_with_real_world_fixture
    kdb = KrakenDbBuilder(db=db)
    
    # run the filter - it should not filter out any sequences because there are no
    # sequences labelled as RSV A in the fixtures
    rsv_min_seq_len= 15000
    kdb._apply_size_filter_to_labelled_sequences(categories=['RSV A'], min_seq_len= rsv_min_seq_len)
    
    update_ids=[3,5,7]
    stmt1="SELECT id, include from sequences WHERE id IN (?,?,?)"
    rows = db._cur.execute(stmt1,update_ids).fetchall()
    excluded_ids = [x['id'] for x in rows if x['include']==0]
    assert not excluded_ids, 'before we apply the RSV filter, none of the three randomly chosen sequences are marked as excluded'
    
    # label some random sequences as "RSV A" (doesn't matter what they really are) 
    # and set their length to slightly below the cutoff so they should be filtered out by the RSB size filter
    stmt2="UPDATE sequences SET category = 'RSV A', seq_length = 14500 WHERE id = ?"
    for id in update_ids:
        db._cur.execute(stmt2,[id])
    db._con.commit()
    
    # apply the filter again - the three sequences should now be filtered out
    kdb._apply_size_filter_to_labelled_sequences(categories=['RSV A'], min_seq_len= rsv_min_seq_len)
    rows = db._cur.execute(stmt1,update_ids).fetchall()
    excluded_ids = [x['id'] for x in rows if x['include']==0]
    assert sorted(excluded_ids) == sorted(update_ids), 'having set the RSV A label and sequence length below cutoff for three randomly chosen sequences, these are now marked as excluded'
    
def test_filter_out_sequences_linked_to_taxonomy_sub_tree(setup_db_with_real_world_fixture):
    """
    TODO: This is highly repetetive with the test test_db.py::test_get_sequence_ids_linked_to_taxon
    because the method tested here is just applying two Db class methods to retrieve, then filter out the
    sequences that are linked to a sub-tree. It would be better to have some fixtures for this that 
    can be shared instead of manipulating the fixtures in the same way here and in the test_db.py method.  
    """
    db = setup_db_with_real_world_fixture
    kdb = KrakenDbBuilder(db=db)
    excluded_seq_stmt = """
        SELECT id
        FROM sequences 
        WHERE id IN (1,2,3,4,6) AND include = 0
    """
    seq_ids = kdb.filter_out_sequences_linked_to_taxonomy_sub_tree(tax_id= 2955291)
    assert not seq_ids, 'before we make changes to the fixtures, no sequences are linked to taxon 2955291 or any node in the sub-tree rooted at this taxon, so none are removed'
    assert not db._cur.execute(excluded_seq_stmt).fetchall(), '... and none of the sequences with ids 1,2,3,4 and 6 are marked include=0'
    
    # For details on these fixture UPDATEs see test_db.py::test_get_sequence_ids_linked_to_taxon
    update_data = [
        [11320,1],
        [114727,2],
        [211044,3],
        [641809,4],
        [335341,6]
    ]
    db._cur.executemany("UPDATE sequences SET tax_id = ? WHERE id = ?", update_data)
    db._con.commit()
    
    # All of the 5 sequences should now be removed by the method because they are now linked to the
    # tree with root tax_id 2955291
    seq_ids = kdb.filter_out_sequences_linked_to_taxonomy_sub_tree(tax_id= 2955291)
    assert sorted(seq_ids) == [1,2,3,4,6], 'after linking 5 sequences to the sub-tree, they are all removed by the method no and their ids are returned'
    rows = db._cur.execute(excluded_seq_stmt).fetchall()
    assert len(rows) == 5, '... and all of the sequences with ids 1,2,3,4 and 6 are marked include=0'
    


def test_create_rsv_taxonomy(setup_db_with_rsv_fixture):
    db = setup_db_with_rsv_fixture
    kdb = KrakenDbBuilder(db=db)
    
    # assert that the data, before we apply the RSV taxonomy modifications, is in a state similar to what 
    # we would have from a default NCBI RefSeq data build 
    stmt1 = """
        SELECT
            fasta_header,
            parent_tax_id,
            parent_tax_names.name AS parent_taxon_name
        FROM sequences
        INNER JOIN taxonomy_nodes ON(taxonomy_nodes.tax_id = sequences.tax_id)
        INNER JOIN taxonomy_names AS parent_tax_names ON(parent_tax_names.tax_id = taxonomy_nodes.parent_tax_id)
        WHERE parent_tax_names.name = ?
    """
    rows = db._cur.execute(stmt1, ['Human respiratory syncytial virus A']).fetchall()
    assert not rows, 'before we start, no sequences are linked to hRSV A'
    
    rows = db._cur.execute(stmt1, ['Human respiratory syncytial virus B']).fetchall()
    assert not rows, 'before we start, no sequences are linked to hRSV B'
    
    # run the RSV creation
    assert kdb.create_rsv_taxonomy(rsv_size_filter=True), 'method returns True'
    
    rows = db._cur.execute(stmt1, ['Human respiratory syncytial virus A']).fetchall()
    assert  len(rows) == 2, 'having built the RSV taxonomy, there are now 2 sequences linked to hRSV A'
    
    rows = db._cur.execute(stmt1, ['Human respiratory syncytial virus B']).fetchall()
    assert len(rows) == 1, 'having built the RSV taxonomy, there is now 1 sequence linked to hRSV B'


def test_filter_out_sequences_linked_to_high_level_rsv_nodes(setup_db_with_rsv_fixture_for_filter):
    db = setup_db_with_rsv_fixture_for_filter
    kdb = KrakenDbBuilder(db=db)
    
    stmt = """
        SELECT
            tax_id,
            include
        FROM sequences 
        WHERE fasta_header IN (
            'Human orthopneumovirus Subgroup B',
            'Human orthopneumovirus Subgroup A',
            'Respiratory syncytial virus, complete genome',
            'some hRSV virus A',
            'known RSV A',
            'known RSV A 2',
            'known RSV B',
            'short RSV A'
        )
    """
    rows = db._cur.execute(stmt).fetchall()
    assert len(rows) == 8, 'there are 7 RSV sequences in total in the fixtures'
    assert len([x for x in rows if x['include']==1]) ==8 ,'before the filter is run, all are marked as included'

    # create the new RSV taxonomy nodes
    # NOTE: this must be done here rather than being in the fixtures as pre-linked sequences
    # because the creation of new taxonomy nodes is needed to set the range of "artificial" 
    # node IDs that need to be skipped
    kdb.create_rsv_taxonomy()

    n_removed = kdb.filter_out_sequences_linked_to_high_level_rsv_nodes()
    
    # only the RefSeq sequence should be removed by the filter
    assert n_removed == 4, 'filter returns correct number of sequences removed'
    
    rows = db._cur.execute(stmt).fetchall()
    assert len(rows) == 8, 'the filter has not changed the number of RSV sequence records'
    assert len([x for x in rows if x['include']==1]) ==4 ,'after the filter, only the 4 pre-labelled RSV sequences linked to hRSV A/B remain included'


def test_deduplicate_sequences(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    kdb = KrakenDbBuilder(db=db)
    
    stmt = """
        SELECT 
            id, 
            include
        FROM sequences
    """
    rows = db._cur.execute(stmt).fetchall()
    assert len([x for x in rows if x['include']==1]) ==35 ,'we are starting with 35 sequence records in the DB that have include=1'
    
    removed_ids = kdb.deduplicate_sequences()
    assert len(removed_ids) == 0, 'before we create duplicates, no sequences are removed as duplicates'
    
    # add two sequences that are duplicates of 2 existing sequences
    update_stmt = """
        INSERT INTO sequences (id,tax_id,fasta_header,dna_sequence,percent_n,seq_length,segment_number,ncbi_acc,flu_name,flu_type,flu_a_h_subtype,flu_a_n_subtype,include,is_flu,category,original_tax_id) VALUES
            (200,NULL,'kraken:taxid|211044|NC_002021.1 Influenza A virus (A/Puerto Rico/8/1934(H1N1)) segment 2, complete sequence','AGCGAAAGCAGGCAAACCATTTGAATGGATGTCAATCCGACCTTACTTTTCTTAAAAGTGCCAGCACAAAATGCTATAAGCACAACTTTCCCTTATACCGGAGACCCTCCTTACAGCCATGGGACAGGAACAGGATACACCATGGATACTGTCAACAGGACACATCAGTACTCAGAAAAGGCAAGATGGACAACAAACACCGAAACTGGAGCACCGCAACTCAACCCGATTGATGGGCCACTGCCAGAAGACAATGAACCAAGTGGTTATGCCCAAACAGATTGTGTATTGGAAGCAATGGCTTTCCTTGAGGAATCCCATCCTGGTATTTTTGAAAACTCGTGTATTGAAACGATGGAGGTTGTTCAGCAAACACGAGTAGACAAGCTGACACAAGGCCGACAGACCTATGACTGGACTTTAAATAGAAACCAGCCTGCTGCAACAGCATTGGCCAACACAATAGAAGTGTTCAGATCAAATGGCCTCACGGCCAATGAGTCTGGAAGGCTCATAGACTTCCTTAAGGATGTAATGGAGTCAATGAAAAAAGAAGAAATGGGGATCACAACTCATTTTCAGAGAAAGAGACGGGTGAGAGACAATATGACTAAGAAAATGATAACACAGAGAACAATAGGTAAAAGGAAACAGAGATTGAACAAAAGGAGTTATCTAATTAGAGCATTGACCCTGAACACAATGACCAAAGATGCTGAGAGAGGGAAGCTAAAACGGAGAGCAATTGCAACCCCAGGGATGCAAATAAGGGGGTTTGTATACTTTGTTGAGACACTGGCAAGGAGTATATGTGAGAAACTTGAACAATCAGGGTTGCCAGTTGGAGGCAATGAGAAGAAAGCAAAGTTGGCAAATGTTGTAAGGAAGATGATGACCAATTCTCAGGACACCGAACTTTCTTTGACCATCACTGGAGATAACACCAAATGGAACGAAAATCAGAATCCTCGGATGTTTTTGGCCATGATCACATATATGACCAGAAATCAGCCCGAATGGTTCAGAAATGTTCTAAGTATTGCTCCAATAATGTTCTCAAACAAAATGGCGAGACTGGGAAAAGGGTATATGTTTGAGAGCAAGAGTATGAAACTTAGAACTCAAATACCTGCAGAAATGCTAGCAAGCATTGATTTGAAATATTTCAATGATTCAACAAGAAAGAAGATTGAAAAAATCCGACCGCTCTTAATAGAGGGGACTGCATCATTGAGCCCTGGAATGATGATGGGCATGTTCAATATGTTAAGCACTGTATTAGGCGTCTCCATCCTGAATCTTGGACAAAAGAGATACACCAAGACTACTTACTGGTGGGATGGTCTTCAATCCTCTGACGATTTTGCTCTGATTGTGAATGCACCCAATCATGAAGGGATTCAAGCCGGAGTCGACAGGTTTTATCGAACCTGTAAGCTACATGGAATCAATATGAGCAAGAAAAAGTCTTACATAAACAGAACAGGTACATTTGAATTCACAAGTTTTTTCTATCGTTATGGGTTTGTTGCCAATTTCAGCATGGAGCTTCCCAGTTTTGGTGTGTCTGGGAGCAACGAGTCAGCGGACATGAGTATTGGAGTTACTGTCATCAAAAACAATATGATAAACAATGATCTTGGTCCAGCAACAGCTCAAATGGCCCTTCAGTTGTTCATCAAAGATTACAGGTACACGTACCGATGCCATAGAGGTGACACACAAATACAAACCCGAAGATCATTTGAAATAAAGAAACTGTGGGAGCAAACCCGTTCCAAAGCTGGACTGCTGGTCTCCGACGGAGGCCCAAATTTATACAACATTAGAAATCTCCACATTCCTGAAGTCTGCCTAAAATGGGAATTGATGGATGAGGATTACCAGGGGCGTTTATGCAACCCACTGAACCCATTTGTCAGCCATAAAGAAATTGAATCAATGAACAATGCAGTGATGATGCCAGCACATGGTCCAGCCAAAAACATGGAGTATGATGCTGTTGCAACAACACACTCCTGGATCCCCAAAAGAAATCGATCCATCTTGAATACAAGTCAAAGAGGAGTACTTGAAGATGAACAAATGTACCAAAGGTGCTGCAATTTATTTGAAAAATTCTTCCCCAGCAGTTCATACAGAAGACCAGTCGGGATATCCAGTATGGTGGAGGCTATGGTTTCCAGAGCCCGAATTGATGCACGGATTGATTTCGAATCTGGAAGGATAAAGAAAGAAGAGTTCACTGAGATCATGAAGATCTGTTCCACCATTGAAGAGCTCAGACGGCAAAAATAGTGAATTTAGCTTGTCCTTCATGAAAAAATGCCTTGTTCCTACT',0,2341,2,'NC_002021.1','A/Puerto Rico/8/1934(H1N1)','A',1,1,1,1,NULL,211044),
            (300,NULL,'kraken:taxid|211044|NC_002022.1 Influenza A virus (A/Puerto Rico/8/1934(H1N1)) segment 3, complete sequence','AGCGAAAGCAGGTACTGATCCAAAATGGAAGATTTTGTGCGACAATGCTTCAATCCGATGATTGTCGAGCTTGCGGAAAAAACAATGAAAGAGTATGGGGAGGACCTGAAAATCGAAACAAACAAATTTGCAGCAATATGCACTCACTTGGAAGTATGCTTCATGTATTCAGATTTCCACTTCATCAATGAGCAAGGCGAGTCAATAATCGTAGAACTTGGTGATCCTAATGCACTTTTGAAGCACAGATTTGAAATAATCGAGGGAAGAGATCGCACAATGGCCTGGACAGTAGTAAACAGTATTTGCAACACTACAGGGGCTGAGAAACCAAAGTTTCTACCAGATTTGTATGATTACAAGGAAAATAGATTCATCGAAATTGGAGTAACAAGGAGAGAAGTTCACATATACTATCTGGAAAAGGCCAATAAAATTAAATCTGAGAAAACACACATCCACATTTTCTCGTTCACTGGGGAAGAAATGGCCACAAAGGCCGACTACACTCTCGATGAAGAAAGCAGGGCTAGGATCAAAACCAGGCTATTCACCATAAGACAAGAAATGGCCAGCAGAGGCCTCTGGGATTCCTTTCGTCAGTCCGAGAGAGGAGAAGAGACAATTGAAGAAAGGTTTGAAATCACAGGAACAATGCGCAAGCTTGCCGACCAAAGTCTCCCGCCGAACTTCTCCAGCCTTGAAAATTTTAGAGCCTATGTGGATGGATTCGAACCGAACGGCTACATTGAGGGCAAGCTGTCTCAAATGTCCAAAGAAGTAAATGCTAGAATTGAACCTTTTTTGAAAACAACACCACGACCACTTAGACTTCCGAATGGGCCTCCCTGTTCTCAGCGGTCCAAATTCCTGCTGATGGATGCCTTAAAATTAAGCATTGAGGACCCAAGTCATGAAGGAGAGGGAATACCGCTATATGATGCAATCAAATGCATGAGAACATTCTTTGGATGGAAGGAACCCAATGTTGTTAAACCACACGAAAAGGGAATAAATCCAAATTATCTTCTGTCATGGAAGCAAGTACTGGCAGAACTGCAGGACATTGAGAATGAGGAGAAAATTCCAAAGACTAAAAATATGxxxxxxxCAAGTCAGCTAAAGTGGGCACTTGGTGAGAACATGGCACCAGAAAAGGTAGACTTTGACGACTGTAAAGATGTAGGTGATTTGAAGCAATATGATAGTGATGAACCAGAATTGAGGTCGCTTGCAAGTTGGATTCAGAATGAGTTCAACAAGGCATGCGAACTGACAGATTCAAGCTGGATAGAGCTTGATGAGATTGGAGAAGATGTGGCTCCAATTGAACACATTGCAAGCATGAGAAGGAATTATTTCACATCAGAGGTGTCTCACTGCAGAGCCACAGAATACATAATGAAGGGGGTGTACATCAATACTGCCTTACTTAATGCATCTTGTGCAGCAATGGATGATTTCCAATTAATTCCAATGATAAGCAAGTGTAGAACTAAGGAGGGAAGGCGAAAGACCAACTTGTATGGTTTCATCATAAAAGGAAGATCCCACTTAAGGAATGACACCGACGTGGTAAACTTTGTGAGCATGGAGTTTTCTCTCACTGACCCAAGACTTGAACCACACAAATGGGAGAAGTACTGTGTTCTTGAGATAGGAGATATGCTTCTAAGAAGTGCCATAGGCCAGGTTTCAAGGCCCATGTTCTTGTATGTGAGGACAAATGGAACCTCAAAAATTAAAATGAAATGGGGAATGGAGATGAGGCGTTGTCTCCTCCAGTCACTTCAACAAATTGAGAGTATGATTGAAGCTGAGTCCTCTGTCAAAGAGAAAGACATGACCAAAGAGTTCTTTGAGAACAAATCAGAAACATGGCCCATTGGAGAGTCTCCCAAAGGAGTGGAGGAAAGTTCCATTGGGAAGGTCTGCAGGACTTTATTAGCAAAGTCGGTATTTAACAGCTTGTATGCATCTCCACAACTAGAAGGATTTTCAGCTGAATCAAGAAAACTGCTTCTTATCGTTCAGGCTCTTAGGGACAATCTGGAACCTGGGACCTTTGATCTTGGGGGGCTATATGAAGCAATTGAGGAGTGCCTAATTAATGATCCCTGGGTTTTGCTTAATGCTTCTTGGTTCAACTCCTTCCTTACACATGCATTGAGTTAGTTGTGGCAGTGCTACTATTTGCTATCCATACTGTCCAAAAAAGTACCTTGTTTCTACT',0,2233,3,'NC_002022.1','A/Puerto Rico/8/1934(H1N1)','A',1,1,1,1,NULL,211044)
    """
    db._cur.executescript(update_stmt)
    db._con.commit()
    
    rows = db._cur.execute(stmt).fetchall()
    assert len([x for x in rows if x['include']==1]) ==37 ,'2 sequences add but still all have include=1'
    
    removed_ids = kdb.deduplicate_sequences()
    assert len(removed_ids) == 2, 'the deduplicate filter now removes 2 sequences'
    removed_ids == [200,300], 'the correct sequences have been reported as removed'
    
    rows = db._cur.execute(stmt).fetchall()
    excluded_ids = [int(x['id']) for x in rows if x['include']==0]
    assert 200 in excluded_ids, 'duplicate seq ID 200 been filtered out'
    assert 300 in excluded_ids, 'duplicate seq ID 200 been filtered out'
        
    included_ids = [int(x['id']) for x in rows if x['include']==1]
    assert 2 in included_ids, 'original seq ID 2 has not been filtered out'
    assert 3 in included_ids, 'original seq ID 3 has not been filtered out'

def test_find_multiref_paths(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    kdb = KrakenDbBuilder(db=db)

    ## mae first call to method
    multiref_paths, seen, multiref_data = kdb.find_multiref_paths()
    ## output should be a dict
    assert isinstance(multiref_paths, list)
    ## output should be empty
    assert len(multiref_paths) == 0

    update_stmt = """
    UPDATE sequences
    SET tax_id = 11250
    WHERE original_tax_id = 11250
    """
    db._cur.execute(update_stmt)
    db._con.commit()

    test_stmt = """
    SELECT tax_id
    FROM sequences
    WHERE original_tax_id = 11250
    """

    ## confirm that the update statement worked
    data = db._cur.execute(test_stmt).fetchall()
    for i in data:
        taxid_found = i[0]
    assert taxid_found == 11250

    ## second call to method
    multiref_paths_now, seen, multiref_data_now = kdb.find_multiref_paths()
    ## output should be a dict
    assert isinstance(multiref_paths_now, list)
    ## added the taxid "11250" to the tax_id field in the sequences table
    ## there are 3 paths in the fixture DB where "11250" is the parent node
    ## output should therefore contain 3 paths
    assert len(multiref_paths_now) == 3
    for path in multiref_paths_now:
        assert 11250 in path

def test_repair_multiref_paths(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    kdb = KrakenDbBuilder(db=db)
    update_stmt_subterminal = """
    UPDATE sequences
    SET tax_id = 11250
    WHERE original_tax_id = 11250
    """

    update_stmt_complex = """
    UPDATE sequences
    SET tax_id = 2955465
    WHERE original_tax_id = 518987
    """
    db._cur.execute(update_stmt_subterminal)
    db._cur.execute(update_stmt_complex)
    db._con.commit()

    multiref_paths, seen, multiref_data = kdb.find_multiref_paths()
    assert len(multiref_paths) == 4

    kdb.repair_multiref_paths(multiref_paths, seen)

    multiref_paths_now, seen, multiref_data_now = kdb.find_multiref_paths()
    assert len(multiref_paths_now) == 0
    
def test_filter_out_sequences_linked_to_high_level_flu_nodes(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    kdb = KrakenDbBuilder(db=db)
    
    # prepare the DB by creating the custom flu taxonomy
    # TODO: should create a DB fixture for this case with custom nodes already created
    kdb.create_segmented_flu_taxonomy_nodes()
    
    n_removed = kdb.filter_out_sequences_linked_to_high_level_flu_nodes()
    assert n_removed == 0, 'no sequences have linkages to any taxa in the fixtures initially, hence nothing is filtered out'
    
    # link:
    # sequence ID1 1 to the custom "Influenza A H1 segment 4" node
    # sequence ID 2 to "H1N1 subtype"
    # -> sequence 1 should remain and sequence 2 should be filtered out after running filter
    inf_a_h1_seg4_tax_id = kdb._db.retrieve_tax_id_by_node_scientific_name("Influenza A H1 segment 4")
    assert inf_a_h1_seg4_tax_id
    h1n1_sub_id = kdb._db.retrieve_tax_id_by_node_scientific_name("H1N1 subtype")
    assert h1n1_sub_id
    
    db._cur.execute("UPDATE sequences SET tax_id = ? WHERE id = ?", [inf_a_h1_seg4_tax_id, 1])
    db._cur.execute("UPDATE sequences SET tax_id = ? WHERE id = ?", [h1n1_sub_id, 2])
    db._con.commit()
    
    n_removed = kdb.filter_out_sequences_linked_to_high_level_flu_nodes()
    assert n_removed == 1, 'having added two sequences to the taxonomy tree (one to a custom node and one bypassing custom taxonomy), one sequence is now filtered out'
    
    row = db._cur.execute("SELECT id, include FROM sequences where id = 1").fetchone()
    assert row['include'] == 1 ,'a sequence that is linked to a new custom flu tax node is still included after the filter'
    
    row = db._cur.execute("SELECT id, include FROM sequences where id = 2").fetchone()
    assert row['include'] == 0 ,'a sequence that is linked to a tax node outside of the custom flu taxonomy (H1N1 subtype) is filtered out after running the filter'
    
def test_filter_max_percent_n(setup_db_with_real_world_fixture):
    db = setup_db_with_real_world_fixture
    kdb = KrakenDbBuilder(db=db)
    
    num_seq_ids_removed = kdb.filter_max_percent_n(10)
    assert num_seq_ids_removed == 0 ,'there are no sequences in the fixtures with Ns, so none are removed'
    
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
    
    select_stmt = "SELECT id, include include FROM sequences where id = 1002"
    row = db._cur.execute(select_stmt).fetchone()
    assert row['include']==1, 'before we apply the filter, seq ID 1002 is included (not filtered out)'
    
    num_seq_ids_removed = kdb.filter_max_percent_n(10)
    
    assert num_seq_ids_removed == 1 ,'there is now one sequence in the DB with >10% N'
    row = db._cur.execute(select_stmt).fetchone()
    assert row['include']==0, 'after we apply the filter, seq ID 1002 is excluded (filtered out)'
    