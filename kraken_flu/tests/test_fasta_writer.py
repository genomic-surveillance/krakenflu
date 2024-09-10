import pytest
import os.path
from importlib_resources import files

from kraken_flu.src.fasta_writer import write_fasta, _add_taxid, _remove_taxid
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
    
def test__add_taxid():
    assert _add_taxid('some random text',1234) == 'kraken:taxid|1234|some random text'
    
def test__remove_taxid():
    assert _remove_taxid('kraken:taxid|1234|some random text') == 'some random text'
    
def test_write_fasta(setup_db_with_real_world_fixture, tmp_path):
    db = setup_db_with_real_world_fixture
    out_file = tmp_path / "out.fna"

    assert not out_file.is_file(), 'before we start, the output file does not exist'
    assert write_fasta(db=db, path= out_file)
    assert out_file.is_file(), 'after running write_fasta, the output file exists now'

    out_file_rows  = [line.rstrip() for line in open( out_file ) ]
    assert len([x for x in out_file_rows if x.startswith(">")]) == 35, '35 records have been written to FASTA file'
    
    header= 'kraken:taxid|335341|NC_007370.1 Influenza A virus (A/New York/392/2004(H3N2)) segment 8, complete sequence' 
    assert [x for x in out_file_rows if x == '>'+header], 'a fasta header form fixtures can be found in the output FASTA file'
    
    # set the tax_id and a mod_fasta_header for the above record and re-export
    stmt="""
    UPDATE sequences
    SET tax_id = ?, mod_fasta_header = ?
    WHERE fasta_header = ?
    """
    new_tax_id= 1234
    mod_fasta_header = 'this is the new header'
    db._cur.execute(stmt,[new_tax_id, mod_fasta_header, header])
    db._con.commit()
    
    assert write_fasta(db=db, path= out_file)
    out_file_rows  = [line.rstrip() for line in open( out_file ) ]
    expected_header = 'kraken:taxid|' + str(new_tax_id) + "|" + mod_fasta_header
    assert '>'+expected_header in out_file_rows, 'after setting tax_id and mod_fasta_header, the output header uses the new information'
    
    # set a sequence as excluded and check it is no longer written to output
    stmt="""
    UPDATE sequences
    SET include = 0
    WHERE fasta_header = ?
    """
    db._cur.execute(stmt,[header])
    db._con.commit()
    
    assert write_fasta(db=db, path= out_file)
    out_file_rows  = [line.rstrip() for line in open( out_file ) ]
    assert not '>'+expected_header in out_file_rows, 'having set a sequence as included=0, it is no longer in the output'
    
    
    