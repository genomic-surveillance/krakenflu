import pytest
import os.path
import re
from importlib_resources import files

from kraken_flu.src.ncbi_influenza_fasta_filter import NcbiInfluenzaFastaFilter

FIXTURE_DIR = files('kraken_flu.tests.fixtures')
NCBI_TEST_FILE = FIXTURE_DIR.joinpath(os.path.join('all_ncbi_flu_download','ncbi_flu.fna'))

def test_init():
    filter = NcbiInfluenzaFastaFilter( fasta_file_path=NCBI_TEST_FILE)
    assert filter
    
def test_matching_names():
    filter = NcbiInfluenzaFastaFilter( fasta_file_path=NCBI_TEST_FILE)
    data = filter.records_by_name
    assert isinstance( data, dict)
    assert len( data.keys()) == 4, 'FASTA fixture files contains data for 20 genomes (distinct virus names)'
    
    name='Influenza B virus (B/Texas/24/2020)'
    assert name in data ,'an expected genome name is found in the data'
    assert isinstance( data[name], dict), '...the dictionary has a 2nd level'
    assert '2' in data[name], '...segment 1 is a 2nd level dict key for this genome'
    assert isinstance( data[name]['2'], dict),'the dict has a 3rd level'
    assert data[name]['2']['fasta_head'] == 'gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds' ,'...header name correct'
    assert  data[name]['2']['seq_len'] == 2367, '...sequence length is correctly extracted'
    
def test_filtered_fasta_headers():
    filter = NcbiInfluenzaFastaFilter( fasta_file_path=NCBI_TEST_FILE)
    data = filter.filtered_fasta_headers
    assert isinstance( data, list)
    assert len( data ) == 2 * 8, 'there are two complete genomes in the fixture file, 2*8 FASTA headers`'
    