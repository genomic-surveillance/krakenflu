import argparse
import os
import shutil
import logging
from kraken_flu.src.kraken_db_builder import KrakenDbBuilder
from kraken_flu.src.ncbi_influenza_fasta_filter import NcbiInfluenzaFastaFilter

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

    subparsers = parser.add_subparsers(help='commands')
    subparsers.required = True
    subparsers.dest = 'command'
    
    # filter command
    filter_parser = subparsers.add_parser(
        'filter', 
        help='filter a large FASTA file of influenza genomes before building a library from it')
    
    # set the function to dispatch to
    filter_parser.set_defaults( func = filter )
    filter_parser.add_argument(
        '--in_fasta',
        required= True,
        metavar= 'FILE',
        type = str,
        help='path to the input FASTA file')
    
    filter_parser.add_argument(
        '--out_fasta',
        required= True,
        metavar= 'FILE',
        type = str,
        help='path to the output (filtered) FASTA file')
    
    filter_parser.add_argument(
        '--discard_duplicates',
        required= False,
        action= 'store_true',
        help='if used, FASTA headers where the name and segment number is not unique are discarded')
    
    build_parser = subparsers.add_parser(
        'build', 
        help='run the main "build" command of the tool, which creates the files for kraken-build')

    build_parser.set_defaults( func = build )

    build_parser.add_argument(
        '--taxonomy_path','-t',
        action = 'store',
        required = True,
        metavar = 'DIR', 
        type = str,
        help='path to the NCBI taxonomy directory that contains files nodes.dmp and names.dmp')

    build_parser.add_argument(
        '--library_path','-l',
        action = 'store',
        required = True,
        metavar = 'DIR', 
        type = str,
        help='path to the NCBI library director that contains the library.fna FASTA file')

    build_parser.add_argument(
        '--acc2tax_path','-a',
        type = str,
        action = 'store',
        default = None,
        metavar = 'FILE',
        required = False,
        help ='path to the NCBI file nucl_gb.accession2taxid IF one was downloaded which is not always the case')

    build_parser.add_argument(
        '--out_dir','-o',
        type = str,
        action = 'store',
        required = True,
        metavar = 'DIR', 
        help = 'path to the output directory which must not exist yet'
    )
    
    return parser
    
def filter(args):
    """
    The filter command
    """
    ncbi_filter = NcbiInfluenzaFastaFilter( fasta_file_path= args.in_fasta, discard_duplicates = args.discard_duplicates )
    ncbi_filter.write_filtered_file( out_path= args.out_fasta )
    
def build(args):
    """
    The build command
    """  
    
    ncbif = KrakenDbBuilder( 
        taxonomy_path= args.taxonomy_path, 
        library_path= args.library_path,
        acc2tax_file_path= args.acc2tax_path)

    ncbif.create_db_ready_dir( args.out_dir )
    
def main():
    
    args = args_parser().parse_args()
    return args.func(args)

if __name__ == "__main__":
    exit(main())