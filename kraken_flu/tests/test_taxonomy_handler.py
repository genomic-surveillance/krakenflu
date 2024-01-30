import pytest
import os.path
import re
from importlib_resources import files

from kraken_flu.src.taxonomy_handler import TaxonomyHandler

FIXTURE_DIR = files('kraken_flu.tests.fixtures')
TAX_DIR = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','taxonomy'))

def test_init():
    th = TaxonomyHandler( taxonomy_path= TAX_DIR)
    assert th
    
    assert os.path.basename( th.names_file_path ) == 'names.dmp', 'found the names file'
    assert os.path.basename( th.nodes_file_path ) == 'nodes.dmp', 'found the nodes file'
    
def test__read_tax_data_file_row():
    th = TaxonomyHandler( taxonomy_path= TAX_DIR)
    row = '11320	|	Influenza A virus	|		|	scientific name	|'
    assert th._read_tax_data_file_row(row) == [ '11320', 'Influenza A virus', '', 'scientific name']
    
def test_names():
    th = TaxonomyHandler( taxonomy_path= TAX_DIR)
    names = th.names
    assert names
    assert isinstance(names, dict)
    assert len(names.keys()) == 21 ,'21 unique tax IDs in the names file'
    
    # use the 'first' (essentially random) key to get a second level element of the data
    names_2nd_lvl = names[list(names.keys())[0]]
    assert isinstance(names_2nd_lvl, list), 'second level is a list'
    
    # 3rd level of the data
    names_3rd_lvl = names[list(names.keys())[0]][0]
    assert isinstance(names_3rd_lvl, dict), 'third level is a dict'
    assert sorted(list(names_3rd_lvl.keys())) == sorted(['name','uname','nclass']), 'third level dict has the expected keys'
    
    assert 11320 in names, '11320 is one of the tax ID keys'
    assert len(names[11320]) == 2, '...there are two names recorded for tax ID 11320'
    assert [x for x in names[11320] if x['name']=='influenza A virus' and x['nclass']=='equivalent name'], 'can find a specific name under tax ID 11320'
    
def test_nodes():
    th = TaxonomyHandler( taxonomy_path= TAX_DIR)
    nodes = th.nodes
    assert nodes
    assert isinstance(nodes, dict)
    assert len(nodes.keys()) == 20, '20 unique tax IDs in the input file (20 rows)'
    
    # use the 'first' (essentially random) key to get a second level element of the data
    nodes_2nd_lvl = nodes[list(nodes.keys())[0]]
    assert isinstance(nodes_2nd_lvl, dict), 'second level is a dict'
    
    assert 11320 in nodes, 'found expected tax ID in nodes'
    assert 'parent_id' in nodes[11320], '...dict has key parent_id'
    assert 'data' in nodes[11320], '...dict has key data'
    assert isinstance(nodes[11320]['data'], list) ,'...data is a list'
    assert len(nodes[11320]['data']) == 11, '...11 items in data list'
    assert nodes[11320]['data'][0] == 'no rank', 'first element of data is the expected value'
    assert nodes[11320]['parent_id'] == 2955291, 'correct parent ID recorded'
    
    assert 11520 in nodes, 'found expected tax ID in nodes'
    assert nodes[11520]['parent_id'] == 2955465, 'correct parent ID recorded'

def test_max_tax_id():
    th = TaxonomyHandler( taxonomy_path= TAX_DIR)
    assert th.max_tax_id() == 3053764, 'found correct max tax ID from the data'

def test_add_taxon():
    th = TaxonomyHandler( taxonomy_path= TAX_DIR)
    
    # create a new taxon as child of 11520 (Influenza B virus)
    max_tax_id_start = th.max_tax_id()
    new_tax_id = max_tax_id_start + 1
    existing_parent_id = 11520
    
    assert new_tax_id not in th.nodes, 'before inserting a new taxon, the tax ID is not present in nodes'
    assert new_tax_id not in th.names, 'before inserting a new taxon, the tax ID is not present in names'
    
    th.add_taxon( tax_id= new_tax_id, parent_tax_id= existing_parent_id, name='a brand new taxon')
    
    assert new_tax_id in th.nodes, 'after inserting a new taxon, the tax ID is present in nodes'
    assert new_tax_id in th.names, 'after inserting a new taxon, the tax ID is present in names'

    assert th.nodes[new_tax_id]['parent_id'] == existing_parent_id, 'the new taxon has correct parent taxon ID'
    assert len(th.names[new_tax_id]) == 1, 'only one name associated with the new taxon ID'
    assert th.names[new_tax_id][0]['name'] == 'a brand new taxon', 'correct name recorded for new taxon`'
    
def test__format_tax_data_file_output_row():
    th = TaxonomyHandler( taxonomy_path= TAX_DIR)
    assert th._format_tax_data_file_output_row( data = [ 1,'test','','abc']) == "1\t|\ttest\t|\t\t|\tabc\t|"
    
def test_write_names_file(tmp_path):
    th = TaxonomyHandler( taxonomy_path= TAX_DIR)
    out_file = tmp_path / "names.dmp"
    
    assert not out_file.is_file(),  'before we start, the output file does not exist'
    th.write_names_file( path= out_file )
    assert out_file.is_file(), '...after calling write_names_file(), the file exists'
    
    original_file_rows = [line.rstrip() for line in open( th.names_file_path) ]
    new_file_rows  = [line.rstrip() for line in open( out_file ) ]
    
    assert len(original_file_rows) == len(new_file_rows), 'input and output file have same number of rows'

    fluB_row_infile = [x for x in original_file_rows if 'Influenza B virus' in x and '11520' in x ]
    assert len(fluB_row_infile) == 1, 'a single row found in input file'
    assert fluB_row_infile[0] == "11520\t|\tInfluenza B virus\t|\t\t|\tscientific name\t|"
    
    fluB_row_outfile = [x for x in new_file_rows if 'Influenza B virus' in x and '11520' in x]
    assert len(fluB_row_outfile) == 1, 'a single row found in output file'
    assert fluB_row_outfile == fluB_row_infile, 'the input and output rows are identical'
    
def test_write_nodes_file(tmp_path):
    th = TaxonomyHandler( taxonomy_path= TAX_DIR)
    out_file = tmp_path / "nodes.dmp"
    
    assert not out_file.is_file(),  'before we start, the output file does not exist'
    th.write_nodes_file( path= out_file )
    assert out_file.is_file(), '...after calling write_nodes_file(), the file exists'
    
    original_file_rows = [line.rstrip() for line in open( th.nodes_file_path) ]
    new_file_rows  = [line.rstrip() for line in open( out_file ) ]
    
    assert len(original_file_rows) == len(new_file_rows), 'input and output file have same number of rows'
    
    row_infile_11250 = [x for x in original_file_rows if re.search('^11250',x) ]
    assert len(row_infile_11250) == 1, 'a single row found in input file'
    assert row_infile_11250[0] == "11250\t|\t3049954\t|\tno rank\t|\t\t|\t9\t|\t1\t|\t1\t|\t1\t|\t0\t|\t1\t|\t1\t|\t0\t|\tcode compliant; specified\t|"
    
    row_outfile_11250 = [x for x in new_file_rows if re.search('^11250',x) ]
    assert len(row_outfile_11250) == 1, 'a single row found in input file'
    assert row_outfile_11250 == row_infile_11250, 'the input and output rows are identical'
    