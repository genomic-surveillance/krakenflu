import pytest
import os.path
from importlib_resources import files

from kraken_flu.src.taxonomy_writer import write_taxonomy, _write_names_file, _write_nodes_file, _format_for_tax_file_output
from kraken_flu.src.db import Db

# set this to True to make the tests print all executed SQL or False to stop that
PRINT_SQL_TRACE=False

@pytest.fixture(scope='function')
def setup_db( tmp_path ):
    db_path = tmp_path / 'kraken_flu.db'
    db = Db(db_path, debug=PRINT_SQL_TRACE)
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
    
def test__format_for_tax_file_output():
    data=[11320,'Influenza A virus',None,'scientific name']
    expected_output = "11320\t|\tInfluenza A virus\t|\t\t|\tscientific name	|"
    assert _format_for_tax_file_output(data) == expected_output, 'a list of column data is formatted to the correct NCBI output format'
        
def test__write_nodes_file(setup_db_with_real_world_fixture, tmp_path):
    db = setup_db_with_real_world_fixture
    out_file = tmp_path / "nodes.dmp"

    assert not out_file.is_file(), 'before we start, the output file does not exist'
    assert _write_nodes_file(db=db, path= out_file)
    
    ############# ADD MORE TESTS 