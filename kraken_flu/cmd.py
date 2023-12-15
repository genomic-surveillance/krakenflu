import argparse
import os
import shutil
import logging
from kraken_flu.src.kraken_db_ncbi_files import KrakenDbNcbiFiles

# get the version number from a file that is created by setuptools_scm 
# when the package is installed. 
try:
    from .version import version as __version__
    from .version import version_tuple
except ImportError:
    __version__ = "unknown version"
    version_tuple = (0, 0, "unknown version")

# NOTE about imports and using this script in develpoment and after production.
# The script uses absolute imports that will only work after being installed by pip.
# To run the script during development, run as follows from top level directory of the repo:
# python -m kraken_db_maker.cmd { OPTIONS }

def args_parser():
    """
    Command line argument parser
    """        
    parser = argparse.ArgumentParser(
        description = "kraken_flu.py: tool to modify KRAKEN2 database build files to handle flu segments as individual genomes") 

    parser.add_argument(
        '--version', '-v',
        action='version', 
        version='kraken_flu ' + __version__ )

    parser.add_argument(
        '--taxonomy_path','-t',
        action = 'store',
        required = True,
        metavar = 'DIR', 
        type = str,
        help='path to the NCBI taxonomy directory that contains files nodes.dmp and names.dmp')

    parser.add_argument(
        '--library_path','-l',
        action = 'store',
        required = True,
        metavar = 'DIR', 
        type = str,
        help='path to the NCBI library director that contains the library.fna FASTA file')

    parser.add_argument(
        '--acc2tax_path','-a',
        type = str,
        action = 'store',
        default = None,
        metavar = 'FILE',
        required = False,
        help ='path to the NCBI file nucl_gb.accession2taxid IF one was dowloaded which is not always the case')

    parser.add_argument(
        '--out_dir','-o',
        type = str,
        action = 'store',
        required = True,
        metavar = 'DIR', 
        help = 'path to the output directory which must not exist yet'
    )
    
    return parser
    
def main():
    
    args = args_parser().parse_args()

    ncbif = KrakenDbNcbiFiles( 
        taxonomy_path= args.taxonomy_path, 
        library_path= args.library_path,
        acc2tax_file_path= args.acc2tax_path)

    path = args.out_dir
    if os.path.exists( path ):
        raise ValueError(f'directory { path } exists already. Will not write into existing directory')
    os.mkdir( path )

    library_path = os.path.join( path , 'library' )
    os.mkdir( library_path )
    taxonomy_path = os.path.join( path, 'taxonomy' )
    os.mkdir( taxonomy_path )

    ncbif.write_modified_fasta_file( os.path.join( library_path, 'library.fna' ))
    ncbif.write_modified_names_files( os.path.join( taxonomy_path, 'names.dmp' ))
    ncbif.write_modified_nodes_files( os.path.join( taxonomy_path, 'nodes.dmp' ))

    if args.acc2tax_path:
        shutil.copyfile( ncbif.acc2tax_file_path, os.path.join( path, 'taxonomy', 'nucl_gb.accession2taxid' ) )

if __name__ == "__main__":
    exit(main())