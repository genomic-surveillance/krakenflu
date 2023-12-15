import pytest
import os.path
import re
from importlib_resources import files

from kraken_db_maker.src.kraken_db_ncbi_files import KrakenDbNcbiFiles

FIXTURE_DIR = files('kraken_db_maker.tests.fixtures')
TAX_DIR = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','taxonomy'))
LIB_DIR = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','library','viral'))

def test_init():
    ncbif = KrakenDbNcbiFiles( taxonomy_path= TAX_DIR, library_path=LIB_DIR)
    assert ncbif
    
    assert ncbif.min_new_tax_id == 3053765, 'minimum new taonx ID correctly identified from names.dmp'
    
    assert os.path.basename( ncbif.names_file_path ) == 'names.dmp', 'found the names file'
    assert os.path.basename( ncbif.nodes_file_path ) == 'nodes.dmp', 'found the nodes file'
    assert os.path.basename( ncbif.fasta_file_path ) == 'library.fna', 'found the fasta file library.fna'
        
def test_flu_genomes_ncbi_to_new_tax_and_parent_ids():
    ncbif = KrakenDbNcbiFiles( taxonomy_path= TAX_DIR, library_path=LIB_DIR)
    
    data = ncbif.flu_genomes_ncbi_to_new_tax_and_parent_ids
    assert data
    assert isinstance( data, dict )
    assert len(data) >0
    
    # NCBI acc IDs for A/Puerto Rico/8/1934(H1N1)
    genome_segments1 = [ 
        'NC_002023.1', 
        'NC_002021.1', 
        'NC_002022.1', 
        'NC_002017.1', 
        'NC_002019.1', 
        'NC_002018.1',
        'NC_002016.1',
        'NC_002020.1' ]
    assert data['NC_002023.1']['new_parent_id'] == 211044, 'the new PARENT tax ID is the orignal taxonomy ID from the FASTA header'
    new_parent_ids = [ data[ x ]['new_parent_id'] for x in genome_segments1 ]
    assert len( new_parent_ids ) == 8 ,'all 8 segments of this flu genome are present in the data dict'
    assert len( set( new_parent_ids ) ) == 1 ,'all segments of this flu genome are assigned the same parent taxon ID'
    new_tax_ids = [ data[ x ]['new_tax_id'] for x in genome_segments1 ]
    assert len( set( new_tax_ids ) ) == 8 ,'all new taxonon IDs for the segments of a single flu genome are different'
    assert min( new_tax_ids ) >= ncbif.min_new_tax_id , 'all new tax IDs are greater than the minimum ID for the new tax IDs range'
    
    assert len( data.keys() ) == 4 * 8 ,'there are 4 complete flu genomes in the data (each has 8 NCBI acc IDs)'

def test_write_modified_fasta_file_kraken_build( tmp_path ):
    ncbif = KrakenDbNcbiFiles( taxonomy_path= TAX_DIR, library_path=LIB_DIR)
    out_file = tmp_path / "fasta_out.fna"
    
    assert not out_file.is_file(), 'before we start, the output file does not exist'
    ncbif.write_modified_fasta_file( out_file )
    assert out_file.is_file(), 'after calling the method, the output file now exists'
    
    with open( out_file ) as f:
        header_rows = [ x.strip() for x in f.readlines() if x[0] == '>' ]

    assert len( header_rows ) == 4 * 8 + 3, 'there are 4 segmented flu genomes and 3 non-segmented genomes (RSV and COVID)'
    
    data = ncbif.flu_genomes_ncbi_to_new_tax_and_parent_ids
    new_tax_id_NC_002023 = data['NC_002023.1']['new_tax_id']
    assert [ x for x in header_rows if re.search( f'kraken:taxid|{new_tax_id_NC_002023}|NC_002023.1', x ) ], 'FASTA header was changed in file and shows new taxon ID'
    
    covid_header = '>kraken:taxid|2697049|NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    assert covid_header in header_rows, 'the COVID FASTA header is unchanged in the file'
    
def test_write_modified_fasta_file_not_kraken_build( tmp_path ):
    # this time use files that are not a pre-built Kraken2 DB file set, so tax IDs
    # are not in the FASTA headers and need to be looked up in acc2tax file from NCBI
    tax_dir = FIXTURE_DIR.joinpath(os.path.join('ncbi_data_not_kraken','taxonomy'))
    lib_dir = FIXTURE_DIR.joinpath(os.path.join('ncbi_data_not_kraken','library','viral'))
    acc2tax_file = FIXTURE_DIR.joinpath(os.path.join('ncbi_data_not_kraken','nucl_gb.accession2taxid'))
    
    ncbif = KrakenDbNcbiFiles( taxonomy_path= tax_dir, library_path= lib_dir, acc2tax_file_path= acc2tax_file )
    out_file = tmp_path / "fasta_out.fna"
    
    assert not out_file.is_file(),  'before we start, the output file does not exist'
    ncbif.write_modified_fasta_file( out_file )
    assert out_file.is_file(), 'after calling the method, the output file now exists'
    
    with open( out_file ) as f:
        header_rows = [ x.strip() for x in f.readlines() if x[0] == '>' ]

    assert len( header_rows ) == 4 * 8 + 3, 'there are 4 segmented flu genomes and 3 non-segmented genomes (RSV and COVID)'
    
    data = ncbif.flu_genomes_ncbi_to_new_tax_and_parent_ids
    new_tax_id_NC_002023 = data['NC_002023.1']['new_tax_id']
    assert [ x for x in header_rows if re.search( f'kraken:taxid|{new_tax_id_NC_002023}|NC_002023.1', x ) ], 'FASTA header was changed in file and shows new taxon ID'
    
    covid_header = '>NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    assert covid_header in header_rows, 'the COVID FASTA header is unchanged in the file'
    
def test_write_modified_names_files( tmp_path ):
    ncbif = KrakenDbNcbiFiles( taxonomy_path= TAX_DIR, library_path=LIB_DIR)
    out_file = tmp_path / "mod_names.dmp"
    
    assert not out_file.is_file(), 'before we start, the output file does not exist'
    ncbif.write_modified_names_files( out_file )
    assert out_file.is_file(), 'after calling the method, the output file now exists'

    with open( ncbif.names_file_path, 'r' ) as fh:
        original_file_rows = fh.readlines()
    
    with open( out_file, 'r' ) as fh:
        mod_file_rows = fh.readlines()
        
    # each flu genome adds 8 segment names to the names file
    # NOTE: the original whole genome names are not removed from the file, they 
    #       are simply not used anymore by any sequences in the FASTA file, which
    #       is ok, since the taxonomy files are expected to contain more nodes than
    #       are being used by sequences in the FASTA
    assert len(mod_file_rows) == len(original_file_rows) + 4 * 8 , '4 flu viruses with 8 segments each were added'

def test_write_modified_nodes_files( tmp_path ):
    ncbif = KrakenDbNcbiFiles( taxonomy_path= TAX_DIR, library_path=LIB_DIR)
    out_file = tmp_path / "mod_nodes.dmp"
    
    assert not out_file.is_file(), 'before we start, the output file does not exist'
    ncbif.write_modified_nodes_files( out_file )
    assert out_file.is_file(), 'after calling the method, the output file now exists'

    with open( ncbif.nodes_file_path, 'r' ) as fh:
        original_file_rows = fh.readlines()
    
    with open( out_file, 'r' ) as fh:
        mod_file_rows = fh.readlines()
        
    # each flu genome adds 8 segments to the nodes file
    # NOTE: see note in test_write_modified_names_files
    assert len(mod_file_rows) == len(original_file_rows) + 4 * 8 , '4 flu viruses with 8 segments each were added'