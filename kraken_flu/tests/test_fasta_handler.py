import pytest
import os.path
import re
from importlib_resources import files

from kraken_flu.src.fasta_handler import FastaHandler, FastaRecord

"""
The fixture data was created by the following commands
python3 kraken_flu/tests/fixtures/kraken_ncbi_data/make_test_fasta.py --in_fasta /lustre/scratch126/gsu/team112/personal/fs5/rvi_dev/krakenDBs/kraken_flu/library_download_all_influenza_ncbi/influenza.fna --out_fasta kraken_flu/tests/fixtures/downsampled_ncbi_all_flu/all_influenza_1pcerent.fasta --downsample 0.01
cat kraken_flu/tests/fixtures/all_ncbi_flu_download/library.fna kraken_flu/tests/fixtures/downsampled_ncbi_all_flu/all_influenza_1pcerent.fasta kraken_flu/tests/fixtures/kraken_ncbi_data/library/viral/library.fna > kraken_flu/tests/fixtures/downsampled_ncbi_all_flu/all_influenza_1pcerent_plus_refseq_plus_completeflu.fasta

This creates a downsampled dataset form the "all NCBI flu" set plus the viral refseq fixture plus three previously extracted full length influenza genomes
to make sure there are some complete genomes in the set (would be unlikely in the downsampled set otherwise)

"""
FIXTURE_DIR = files('kraken_flu.tests.fixtures')
FASTA_FILE = FIXTURE_DIR.joinpath(os.path.join('downsampled_ncbi_all_flu','all_influenza_1pcerent_plus_refseq_plus_completeflu.fasta'))

def test_init():
    fh = FastaHandler( fasta_file_path=FASTA_FILE)
    assert fh

def test__parse_header():
    fp = FastaHandler( fasta_file_path=FASTA_FILE)
    
    header = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    flu_type, ncbi_acc, kraken_taxid, is_flu, is_fluA, flu_isolate_name, flu_segment_number, h_subtype, n_subtype = fp._parse_header(header)

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
    flu_type, ncbi_acc, kraken_taxid, is_flu, is_fluA, flu_isolate_name, flu_segment_number, h_subtype, n_subtype = fp._parse_header(header)
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
    flu_type, ncbi_acc, kraken_taxid, is_flu, is_fluA, flu_isolate_name, flu_segment_number, h_subtype, n_subtype = fp._parse_header(header)
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
    flu_type, ncbi_acc, kraken_taxid, is_flu, is_fluA, flu_isolate_name, flu_segment_number, h_subtype, n_subtype = fp._parse_header(header)
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

    
def test_data():
    fh = FastaHandler( fasta_file_path=FASTA_FILE)
    data = fh.data
    assert data
    assert isinstance(data, list), 'data is a list'
    assert len(data)>0, 'there are elements in the list'
    assert isinstance(data[0], FastaRecord), 'the first element is a dict'
    assert len(data) == 5505 ,'number of extracted FASTA header data matches number of sequences in the file'
    
    data_CY030663 = [ x for x in data if x.ncbi_acc =='CY030663']
    assert len(data_CY030663) == 1, 'one record with NCBI acc CY030663 in the data'
    assert data_CY030663[0].orig_head == 'gi|169731751|gb|CY030663|Influenza B virus (B/Tennessee/UR06-0431/2007) segment 1, complete sequence'
    assert data_CY030663[0].is_flu == True
    assert data_CY030663[0].flu_seg_num == 1
    assert data_CY030663[0].taxid is None
    assert data_CY030663[0].flu_name == 'B/Tennessee/UR06-0431/2007'
    assert data_CY030663[0].include_in_output == True ,'filter set to True'
    
    data_NC_002021 = [ x for x in data if x.ncbi_acc =='NC_002021.1']
    assert len(data_NC_002021) == 1, 'one record with NCBI acc NC_002021.1 in the data'
    assert data_NC_002021[0].orig_head == 'kraken:taxid|211044|NC_002021.1 Influenza A virus (A/Puerto Rico/8/1934(H1N1)) segment 2, complete sequence'
    assert data_NC_002021[0].is_flu == True
    assert data_NC_002021[0].flu_seg_num == 2
    assert data_NC_002021[0].taxid == 211044
    assert data_NC_002021[0].flu_name == 'A/Puerto Rico/8/1934(H1N1)'
    
    # this is an unusual case from RefSeq: the influenza B segment 1 is
    # recorded with a non-standard header:
    # >kraken:taxid|518987|NC_002204.1 Influenza B virus RNA 1, complete sequence
    # the isolate name should be derived from other headers by looking up the taxid
    data_NC_002204 = [ x for x in data if x.ncbi_acc == 'NC_002204.1']
    assert len(data_NC_002204) == 1, 'one record with NCBI acc NC_002204.1 in the data'
    assert data_NC_002204[0].orig_head == 'kraken:taxid|518987|NC_002204.1 Influenza B virus RNA 1, complete sequence'
    assert data_NC_002204[0].is_flu == True
    assert data_NC_002204[0].flu_seg_num == 1 ,'correctly parsed segment number given in non-standard form'
    assert data_NC_002204[0].taxid == 518987
    assert data_NC_002204[0].flu_name == 'B/Lee/1940', 'correctly inferred isolate name by looking up taxid in other records (isolate name missing from this header)'

    # make sure that non-flu records are also included in the data
    data_covid = [ x for x in data if x.ncbi_acc == 'NC_045512.2']
    assert len(data_covid) == 1, 'one record with NCBI acc for the covid ref genome is in the data'
    assert data_covid[0].orig_head == 'kraken:taxid|2697049|NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    assert data_covid[0].is_flu == False, 'this is covid, not flu'
    assert data_covid[0].flu_seg_num is None ,"no flu segment number because it isn't flu'"
    assert data_covid[0].taxid == 2697049



def test_remove_incomplete_flu():
    # using a different (smaller) dataset for testing this method
    fh = FastaHandler( fasta_file_path= FIXTURE_DIR.joinpath(os.path.join('all_ncbi_flu_download','library.fna') ) )
    
    assert len( fh.data ) == 20, 'there are 20 records in the FASTA file in total'
    assert len( [ x for x in fh.data if not x.include_in_output ] ) == 0 , 'there are no filtered records before applying the filter'
    
    assert fh.remove_incomplete_flu() == True , 'method returns True'
    
    assert len( [ x for x in fh.data if not x.include_in_output ] ) == 4 ,'after applying the filter, there are 4 filtered records (incomplete flu)'
    assert len( [ x for x in fh.data if x.include_in_output ] ) == 16 ,'after applying the filter, there are 16 complete flu records (not filtered)'

    # check with the bigger dataset
    fh = FastaHandler( fasta_file_path=FASTA_FILE)
    assert len( fh.data ) == 5505, 'there are 5505 records in the FASTA file in total'
    assert len( [ x for x in fh.data if not x.include_in_output ] ) == 0 , 'there are no filtered records before applying the filter'
    
    assert fh.remove_incomplete_flu() == True , 'method returns True'
    assert len( fh.data ) == 5505, 'there are still 5505 records in the FASTA file after filtering (the filter does not remove anything, just marks records)'


def test_n_seq_total():
    fh = FastaHandler( fasta_file_path=FASTA_FILE)
    assert fh.n_seq_total() == 5505 ,'returns correct number of total sequence records'

def test_n_seq_flu():
    fh = FastaHandler( fasta_file_path=FASTA_FILE)
    assert fh.n_seq_flu() == 5502, 'correct number of flu sequences'
    
def test_n_seq_filtered():
    fh = FastaHandler( fasta_file_path= FIXTURE_DIR.joinpath(os.path.join('all_ncbi_flu_download','library.fna') ) )
    assert fh.n_seq_filtered() == 0 ,'before applying filter, n_seq_filtered returns 0'
    assert fh.remove_incomplete_flu(), 'apply flu complete filter'
    assert fh.n_seq_filtered() == 4 ,'correct number of filtered records after applying the filter'


def test_write_fasta( tmp_path ):
    fh = FastaHandler( fasta_file_path=FASTA_FILE)
    out_file = tmp_path / "out.fna"
    
    assert not out_file.is_file(), 'before we start, the output file does not exist'
    fh.write_fasta( out_file )
    assert out_file.is_file(), 'after calling the method, the output file now exists'
    
    out_file_rows  = [line.rstrip() for line in open( out_file ) ]
    assert len([ x for x in out_file_rows if '>' in x ]) == 5505, 'there are 5505 sequence in the output file (=number of seq in input file)'
    
    expected_header = '>kraken:taxid|2697049|NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    assert expected_header in out_file_rows, 'header for SARS-CoV2 should be unmodified in the output file'

    expected_header_regex = re.compile(r'>kraken:taxid\|[0-9]+\|NC_026435.1 Influenza A \(A/California/07/2009\(H1N1\)\) segment 2')
    matching_rows = [ x for x in out_file_rows if expected_header_regex.search( x )]
    print( matching_rows)
    assert len(matching_rows) == 1, 'there is a row in the output that is the expected reformatted header for (A/California/07/2009(H1N1) segment 2'