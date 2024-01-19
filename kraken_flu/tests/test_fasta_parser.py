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
    
    