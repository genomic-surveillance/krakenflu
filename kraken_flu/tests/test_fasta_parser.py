import pytest
import os.path
from importlib_resources import files

from kraken_flu.src.fasta_parser import FastaParser

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
    fp = FastaParser( fasta_file_path=FASTA_FILE)
    assert fp
    
def test_NCBI_ACC_REGEX():
    str = 'gi|59896308|gb|CY000029|Influenza A virus (A/New York/32/2003(H3N2)) segment 8, complete sequence'
    assert FastaParser.NCBI_ACC_REGEX.search(str), 'regex matches expected pattern'
    match = FastaParser.NCBI_ACC_REGEX.search(str)
    assert match.group(0) == 'gb|CY000029', 'regex extracts the correct pattern'
    
    str = '>kraken:taxid|2697049|NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    assert FastaParser.NCBI_ACC_REGEX.search(str), 'regex matches expected pattern'
    match = FastaParser.NCBI_ACC_REGEX.search(str)
    assert match.group(0) == 'NC_045512.2', 'regex extracts the correct pattern'
    
def test_FLU_REGEX():
    str = 'gi|59896308|gb|CY000029|Influenza A virus (A/New York/32/2003(H3N2)) segment 8, complete sequence'
    assert FastaParser.FLU_REGEX.search(str), 'regex matches expected pattern'
    
    str = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    assert FastaParser.FLU_REGEX.search(str), 'regex matches expected pattern'
    
def test_FLU_ISOLATE_NAME_REGEX():
    str = 'gi|59896308|gb|CY000029|Influenza A virus (A/New York/32/2003(H3N2)) segment 8, complete sequence'
    assert FastaParser.FLU_ISOLATE_NAME_REGEX.search(str), 'regex matches expected pattern'
    match = FastaParser.FLU_ISOLATE_NAME_REGEX.search(str)
    assert match.group(1) == 'A/New York/32/2003(H3N2)', 'regex extracts the correct pattern'
    
    str = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    assert FastaParser.FLU_ISOLATE_NAME_REGEX.search(str), 'regex matches expected pattern'
    match = FastaParser.FLU_ISOLATE_NAME_REGEX.search(str)
    assert match.group(1) == 'B/Texas/24/2020', 'regex extracts the correct pattern'
    
    str = '>kraken:taxid|2697049|NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    assert FastaParser.FLU_ISOLATE_NAME_REGEX.search(str) is None, 'does not match'
    
def test_FLU_SEG_NUM_REGEX():
    str = 'gi|59896308|gb|CY000029|Influenza A virus (A/New York/32/2003(H3N2)) segment 8, complete sequence'
    assert FastaParser.FLU_SEG_NUM_REGEX.search(str), 'regex matches expected pattern'
    match = FastaParser.FLU_SEG_NUM_REGEX.search(str)
    assert match.group(1) == '8', 'regex extracts the correct pattern'
    
    str = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    assert FastaParser.FLU_SEG_NUM_REGEX.search(str), 'regex matches expected pattern'
    match = FastaParser.FLU_SEG_NUM_REGEX.search(str)
    assert match.group(1) == '2', 'regex extracts the correct pattern'
    
    str = '>kraken:taxid|2697049|NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    assert FastaParser.FLU_SEG_NUM_REGEX.search(str) is None, 'does not match'
    
    str = 'segment 1 of something that is not flu'
    assert FastaParser.FLU_SEG_NUM_REGEX.search(str) is None, 'does not match'
    
def test_KRAKEN_TAX_ID_REGEX():
    str = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    assert FastaParser.KRAKEN_TAX_ID_REGEX.search(str) is None, 'does not match'
    
    str = '>kraken:taxid|2697049|NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    assert FastaParser.KRAKEN_TAX_ID_REGEX.search(str), 'regex matches expected pattern'
    match = FastaParser.KRAKEN_TAX_ID_REGEX.search(str)
    assert match.group(1) == '2697049', 'regex extracts the correct pattern'

def test__parse_header():
    fp = FastaParser( fasta_file_path=FASTA_FILE)
    
    header = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    ncbi_acc, kraken_taxid, is_flu, flu_isolate_name, flu_segment_number = fp._parse_header(header)
    assert ncbi_acc == 'MT375832', 'correctly parsed NCBI accession from gb|xxxxx format'
    assert kraken_taxid is None ,'no tax ID in this header'
    assert flu_isolate_name == 'B/Texas/24/2020' ,'flu isolate name parsed correctly'
    assert flu_segment_number == 2 ,'correctly parsed flu segment number'
    
    # this time without the leading ">", which is not required for the regexes to work
    header = 'kraken:taxid|641809|NC_026438.1 Influenza A virus (A/California/07/2009(H1N1)) segment 1 polymerase PB2 (PB2) gene, complete cds'
    ncbi_acc, kraken_taxid, is_flu, flu_isolate_name, flu_segment_number = fp._parse_header(header)
    assert ncbi_acc == 'NC_026438.1', 'correctly parsed NCBI accession from non-gb|xxxxx format'
    assert kraken_taxid == 641809 ,'correctly parsed a kraken taxid from this header'
    assert flu_isolate_name == 'A/California/07/2009(H1N1)' ,'flu isolate name parsed correctly'
    assert flu_segment_number == 1 ,'correctly parsed flu segment number'
    
    