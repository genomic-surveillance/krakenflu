import pytest
import os.path
import re
from importlib_resources import files

from kraken_flu.src.kraken_db_ncbi_files import KrakenDbNcbiFiles

FIXTURE_DIR = files('kraken_flu.tests.fixtures')
TAX_DIR = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','taxonomy'))
LIB_DIR = FIXTURE_DIR.joinpath(os.path.join('kraken_ncbi_data','library','viral'))

def test_init():
    ncbif = KrakenDbNcbiFiles( taxonomy_path= TAX_DIR, library_path=LIB_DIR)
    assert ncbif
    
    assert ncbif.min_new_tax_id == 3053765, 'minimum new taonx ID correctly identified from names.dmp'
    
    assert os.path.basename( ncbif.names_file_path ) == 'names.dmp', 'found the names file'
    assert os.path.basename( ncbif.nodes_file_path ) == 'nodes.dmp', 'found the nodes file'
    assert os.path.basename( ncbif.fasta_file_path ) == 'library.fna', 'found the fasta file library.fna'
    assert ncbif.prelim_map_file_path is not None
    assert os.path.basename( ncbif.prelim_map_file_path ) == 'prelim_map.txt'
        
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
    assert [ x for x in header_rows if re.search( rf'>kraken:taxid\|{new_tax_id_NC_002023}\|NC_002023\.1', x ) ], 'FASTA header was changed in file and shows new taxon ID'
    
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
    
def test_flu_genomes_ncbi_to_new_tax_and_parent_ids_all_ncbi_download():
    """
    tests lu_genomes_ncbi_to_new_tax_and_parent_ids against an extract from the download of all
    influenza genomes from NCBI (https://ftp.ncbi.nih.gov/genomes/INFLUENZA/)
    The main difference is that the FASTA headers do not have NC_xxxxx genbank IDs so the IDs have 
    to be extracted from 
    
    """
    lib_dir = FIXTURE_DIR.joinpath(os.path.join('all_ncbi_flu_download' ))
    acc2taxid_file = FIXTURE_DIR.joinpath(os.path.join('all_ncbi_flu_download', 'nucl_gb.accession2taxid' ))
    ncbif = KrakenDbNcbiFiles( taxonomy_path= TAX_DIR, library_path=lib_dir, acc2tax_file_path= acc2taxid_file)
    
    data = ncbif.flu_genomes_ncbi_to_new_tax_and_parent_ids
    assert data
    assert 'MT375832' in data, "a 'gb|xxxxx' formatted Genbank ID was parsed correctly"  
    
def test_write_modified_names_files( tmp_path ):
    ncbif = KrakenDbNcbiFiles( taxonomy_path= TAX_DIR, library_path=LIB_DIR)
    out_file = tmp_path / "mod_names.dmp"
    
    assert not out_file.is_file(), 'before we start, the output file does not exist'
    ncbif.write_modified_names_files( out_file )
    assert out_file.is_file(), 'after calling the method, the output file now exists'

    original_file_rows = [line.rstrip() for line in open( ncbif.names_file_path ) ]
    mod_file_rows  = [line.rstrip() for line in open( out_file ) ]
        
    # each flu genome adds 8 segment names to the names file
    # NOTE: the original whole genome names are not removed from the file, they 
    #       are simply not used anymore by any sequences in the FASTA file.
    #       This is intended behaviour. It does not matter to have more entries in the
    #       taxonomy files than are used in the FASTA because this often happens in kraken
    #       DB creation (the taxonomy files are for all of NCBI)
    assert len(mod_file_rows) == len(original_file_rows) + 4 * 8 , '4 flu viruses with 8 segments each were added'
    known_name_pattern = r'\t\|\tNC_002023\.1 Influenza A virus \(A\/Puerto Rico\/8\/1934\(H1N1\)\) segment 1, complete sequence'
    assert [ r for r in mod_file_rows if re.search( known_name_pattern, r)] , 'we find a row with the exact name of one segment that we expect'


def test_write_modified_nodes_files( tmp_path ):
    ncbif = KrakenDbNcbiFiles( taxonomy_path= TAX_DIR, library_path=LIB_DIR)
    out_file = tmp_path / "mod_nodes.dmp"
    
    assert not out_file.is_file(), 'before we start, the output file does not exist'
    ncbif.write_modified_nodes_files( out_file )
    assert out_file.is_file(), 'after calling the method, the output file now exists'

    original_file_rows = [line.rstrip() for line in open( ncbif.nodes_file_path ) ]
    mod_file_rows  = [line.rstrip() for line in open( out_file ) ]
        
    # each flu genome adds 8 segments to the nodes file
    # NOTE: see note in test_write_modified_names_files
    assert len(mod_file_rows) == len(original_file_rows) + 4 * 8 , '4 flu viruses with 8 segments each were added'
    first_tax_id_pattern = rf'^{ncbif.min_new_tax_id}\t\|\t[0-9]+'
    assert [r for r in mod_file_rows if re.search( first_tax_id_pattern, r)] ,'we find a node with the expected pattern of the first new tax ID in col 1'
    
def test_write_prelim_map_file( tmp_path ):
    ncbif = KrakenDbNcbiFiles( taxonomy_path= TAX_DIR, library_path=LIB_DIR)
    out_file = tmp_path / "prelim_map.txt"
    
    assert not out_file.is_file(), 'before we start, the output file does not exist'
    ncbif.write_prelim_map_file( out_file )
    assert out_file.is_file(), 'after calling the method, the output file now exists'

    original_file_rows = [line.rstrip() for line in open( ncbif.prelim_map_file_path ) ]
    mod_file_rows  = [line.rstrip() for line in open( out_file ) ]
        
    assert len(mod_file_rows) == len(original_file_rows) , 'the prelim_map.txt file has the same number of entries before and after modification'
    
    # the modified file should have the new tax ID twice: once in col 2 as a kraken taxid tag
    # and once on its own in col3 (tab delimited)
    data = ncbif.flu_genomes_ncbi_to_new_tax_and_parent_ids
    new_tax_id_NC_002023 = data['NC_002023.1']['new_tax_id']
    prelim_map_pattern = rf'\tkraken:taxid|{new_tax_id_NC_002023}|NC_002023.1\t{new_tax_id_NC_002023}'
    assert [ x for x in mod_file_rows if re.search( prelim_map_pattern, x ) ], 'the tax ID in the prelim mapping file was changed correctly for this influenza NCBI ID'
    

def test_create_db_ready_dir( tmp_path ):
    out_dir = tmp_path / 'out_dir'
    ncbif = KrakenDbNcbiFiles( taxonomy_path= TAX_DIR, library_path=LIB_DIR)
    
    assert not out_dir.is_dir() , 'before we begin, the output dir does not exist'
    ncbif.create_db_ready_dir( out_dir )
    
    assert out_dir.is_dir(), 'after the command, the out dir exists'
    assert os.path.exists( os.path.join( out_dir, 'library' ) ), 'the library subdir created'
    assert os.path.exists( os.path.join( out_dir, 'taxonomy' ) ), 'the taxonomy subdir created'
    assert os.path.exists( os.path.join( out_dir, 'library', 'library.fna' ) ), 'the FASTA file was created'
    assert os.path.exists( os.path.join( out_dir, 'library', 'prelim_map.txt' ) ), 'the preliminary acc2tax ID file was created'
    assert os.path.exists( os.path.join( out_dir, 'taxonomy', 'names.dmp' ) ), 'the names file was created'
    assert os.path.exists( os.path.join( out_dir, 'taxonomy', 'nodes.dmp' ) ), 'the nodes file was created'

