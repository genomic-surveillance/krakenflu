import os.path
from importlib_resources import files

from kraken_flu.src.kraken_db_builder import KrakenDbBuilder
from kraken_flu.src.db import Db

FIXTURE_DIR = files('kraken_flu.tests.fixtures')
TAX_DIR = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','taxonomy'))

SMALL_VIRUS_FILE = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','library','viral','library.fna'))
OTHER_SMALL_VIRUS_FILE = FIXTURE_DIR.joinpath(os.path.join('ncbi_data_not_kraken','library','viral','library.fna'))

def test_ini_no_patht():
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