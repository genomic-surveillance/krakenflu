import pytest
import os.path
from importlib_resources import files

from kraken_flu.src.fasta_loader import load_fasta, _parse_header, _get_num_records, _calculate_percent_n
from kraken_flu.src.db import Db

FIXTURE_DIR = files('kraken_flu.tests.fixtures')
SMALL_VIRUS_FILE = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','library','viral','library.fna'))

@pytest.fixture(scope='function')
def setup_db( tmp_path ):
    """ initialise the DB """
    db_path = tmp_path / 'kraken_flu.db'
    db = Db(db_path)
    yield db
    
def test__parse_header():
    
    header = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    flu_type, ncbi_acc, kraken_taxid, is_flu, is_fluA, flu_isolate_name, flu_segment_number, h_subtype, n_subtype = _parse_header(header)

    assert ncbi_acc == 'MT375832', 'correctly parsed NCBI accession from gb|xxxxx format'
    assert kraken_taxid is None ,'no tax ID in this header'
    assert is_flu == True ,'this is a flu genome'
    assert flu_type == 'B', 'flu B'
    assert is_fluA == False , 'not a flu type A'
    assert h_subtype == None, 'no flu A so no H subtype'
    assert n_subtype == None , 'no flu A so no N subtype'
    assert flu_isolate_name == 'B/Texas/24/2020' ,'flu isolate name parsed correctly'
    assert flu_segment_number == 2 ,'correctly parsed flu segment number'
    
    # this time without the leading ">", which is not required for the regexes to work
    header = 'kraken:taxid|641809|NC_026438.1 Influenza A virus (A/California/07/2009(H1N1)) segment 1 polymerase PB2 (PB2) gene, complete cds'
    flu_type, ncbi_acc, kraken_taxid, is_flu, is_fluA, flu_isolate_name, flu_segment_number, h_subtype, n_subtype = _parse_header(header)
    assert ncbi_acc == 'NC_026438.1', 'correctly parsed NCBI accession from non-gb|xxxxx format'
    assert kraken_taxid == 641809 ,'correctly parsed a kraken taxid from this header'
    assert is_flu == True ,'this is a flu genome'
    assert flu_type == 'A', 'flu A'
    assert is_fluA == True , 'a flu type A'
    assert h_subtype == 'H1', 'H1 subtype'
    assert n_subtype == 'N1' , 'N1 subtype'
    assert flu_isolate_name == 'A/California/07/2009(H1N1)' ,'flu isolate name parsed correctly'
    assert flu_segment_number == 1 ,'correctly parsed flu segment number'
    
    
    header = '>gi|169731751|gb|CY030663|Influenza B virus (B/Tennessee/UR06-0431/2007) segment 1, complete sequence'
    flu_type, ncbi_acc, kraken_taxid, is_flu, is_fluA, flu_isolate_name, flu_segment_number, h_subtype, n_subtype = _parse_header(header)
    assert ncbi_acc == 'CY030663', 'correctly parsed NCBI accession from gb|xxxxx format'
    assert kraken_taxid is None ,'no tax ID in this header'
    assert is_flu == True ,'this is a flu genome'
    assert flu_type == 'B' ,'flu B'
    assert is_fluA == False , 'not a flu type A'
    assert h_subtype == None, 'no flu A so no H subtype'
    assert n_subtype == None , 'no flu A so no H subtype'
    assert flu_isolate_name == 'B/Tennessee/UR06-0431/2007' ,'flu isolate name parsed correctly'
    assert flu_segment_number == 1 ,'correctly parsed flu segment number'
    
    # a real FASTA header for a Flu B RefSeq - uses the keyword RNA instead of segment
    header = 'kraken:taxid|518987|NC_002204.1 Influenza B virus RNA 1, complete sequence'
    flu_type, ncbi_acc, kraken_taxid, is_flu, is_fluA, flu_isolate_name, flu_segment_number, h_subtype, n_subtype = _parse_header(header)
    assert is_flu == True , 'identified as flu'
    assert ncbi_acc == 'NC_002204.1', 'correctly parsed NCBI accession from RefSeq ID format'
    assert kraken_taxid == 518987 ,'correctly parsed existing kraken tax ID'
    assert is_flu == True ,'this is a flu genome'
    assert is_fluA == False , 'not a flu type A'
    assert flu_type == 'B' , 'flu B'
    assert h_subtype == None, 'no flu A so no H subtype'
    assert n_subtype == None , 'no flu A so no H subtype'
    assert flu_isolate_name is None ,'there is no flu isolate name given in this header'
    assert flu_segment_number == 1 ,'correctly parsed flu segment number'

    
def test_load_fasta( setup_db ):
    db = setup_db
    rows = db._cur.execute("SELECT * FROM sequences").fetchall()
    assert len(rows) == 0, 'no sequence data in DB before we start uploading'
    
    # trigger the DB loading
    assert load_fasta(db, SMALL_VIRUS_FILE), 'loading into DB returns True'
    
    rows = db._cur.execute("SELECT * FROM sequences").fetchall()
    assert len(rows) == 35, '35 rows of sequence data in DB after loadingg'
    
    rows = db._cur.execute("SELECT * FROM sequences WHERE segment_number == 1").fetchall()
    assert len(rows) == 4, '4 sequence records with segment number == 1 (one has RNA 1 instead of segment 1 in the header)'
    
    rows = db._cur.execute("SELECT * FROM sequences WHERE ncbi_acc == 'NC_002211.1'").fetchall()
    assert len(rows) == 1, '1 sequence records with NCBI accession number NC_002211.1'
    
    # retrieve this record and test all inserted fields are as expected
    # >kraken:taxid|335341|NC_007373.1 Influenza A virus (A/New York/392/2004(H3N2)) segment 1, complete sequence
    rows = (
        db._cur.execute(
        "SELECT * FROM sequences WHERE fasta_header == 'kraken:taxid|335341|NC_007373.1 Influenza A virus (A/New York/392/2004(H3N2)) segment 1, complete sequence'")
    .fetchall() )
    assert len(rows) == 1, 'single row returned'
    row = rows[0]
    assert row['ncbi_acc'] == 'NC_007373.1','retrieved row has the correct NCBI accession'
    assert row['dna_sequence'].startswith('AGCAAAAGCA'), '... DNA seq starts is correct'
    assert row['dna_sequence'].endswith('CTTGTTTCTACT'),'... DNA seq ends is correct'
    assert row['seq_length'] == 2341,'... DNA seq length is correct'
    assert row['category'] == None, '... no category is set'
    assert row['flu_type'] == 'A', '... flu type is correct'
    assert row['original_tax_id'] == 335341,'... original kraken taxid is correct'
    assert row['is_flu'], '... is_flu correctly set to True'
    assert row['flu_name'] == 'A/New York/392/2004(H3N2)', '... correct flu isolate name'
    assert row['segment_number'] == 1, '... correct segment number'
    assert row['flu_a_h_subtype'] == 3, '... correct flu H subtype'
    assert row['flu_a_n_subtype'] == 2, '... correct flu N subtype'
    
def test__get_num_records():
    assert _get_num_records(SMALL_VIRUS_FILE) == 35, '35 FASTA records in the file'
    
def test__calculate_percent_n():
    assert isinstance(_calculate_percent_n('nnat'), float)
    assert _calculate_percent_n('nnat') == 50, '50% n detected from lower case sequence string with length calculated on the fly'
    assert _calculate_percent_n('NNaT') == 50, 'upper and lower case letters both work'
    assert _calculate_percent_n('NNATNGCN') == 50, 'Ns can be anywhere in the sequence'
