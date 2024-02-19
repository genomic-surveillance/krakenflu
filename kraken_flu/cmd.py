import argparse
import os
import shutil
import logging
from kraken_flu.src.kraken_db_builder import KrakenDbBuilder

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
# python -m kraken_flu.cmd { OPTIONS }

def args_parser():
    """
    Command line argument parser
    """        
    parser = argparse.ArgumentParser(
        description = "kraken_flu.py: tool to modify KRAKEN2 database build files to handle flu segments as individual genomes") 

    parser.add_argument(
        '-v', '--version', 
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
        '--fasta_path','-l',
        action = 'store',
        required = True,
        metavar = 'DIR', 
        type = str,
        help='path to the sequences FASTA path')

    parser.add_argument(
        '--out_dir','-o',
        type = str,
        action = 'store',
        required = True,
        metavar = 'DIR', 
        help = 'path to the output directory which must not exist yet'
    )
    
    parser.add_argument(
        '--filter','-f',
        action = 'store_true',
        help = 'when used, the FASTA file will be filtered to remove incomplete influenza genomes'        
    )
    
    return parser

def main():
    
    args = args_parser().parse_args()
    kdb = KrakenDbBuilder(
        taxonomy_path= args.taxonomy_path, 
        fasta_file_path= args.fasta_path,
    )    
    kdb.create_db_ready_dir( 
        path = args.out_dir,
        force= True, 
        fasta_file_name= 'library.fna', 
        filter_incomplete_flu= args.filter)

if __name__ == "__main__":
    exit(main())