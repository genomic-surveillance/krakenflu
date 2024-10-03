import argparse
import os
import tempfile
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
        '--no-acc2taxid',
        action = 'store_true',
        help = 'ignore the NCBI acc2taxid file, even if present in the taxonomy path'        
    )

    parser.add_argument(
        '--fasta_path','-f',
        action = 'store',
        required = True,
        metavar = 'FILE', 
        type = str,
        nargs='+',
        help='one or more filepath for sequence FASTA file(s)')

    parser.add_argument(
        '--out_dir','-o',
        type = str,
        action = 'store',
        required = True,
        metavar = 'DIR', 
        help = 'path to the output directory which must not exist yet'
    )
    
    parser.add_argument(
        '--db_file','-d',
        type = str,
        action = 'store',
        required = False,
        metavar = 'FILE', 
        help = 'File path for the sqlite backend DB file created by this tool. If not provided, a tmp file will be used.'
    )

    parser.add_argument(
        '--keep_db_file',
        action = 'store_true',
        help = 'when used, the sqlite backend DB file is not deleted at the end of the process and can be further investigated using sqlite'        
    )
    
    parser.add_argument(
        '--filter_flu',
        action = 'store_true',
        help = 'apply the influenza filters: remove incomplete genomes (also see filter_except) and sequences with non-standard isolate names'        
    )
    
    parser.add_argument(
        '--filter_except',
        action= 'append',
        metavar= 'PATTERN',
        required= False,
        type= str,
        help = 'one or more strings/patterns that are used to exclude genomes from the Influenza "complete genome" filter (if used)'
    )

    parser.add_argument(
        '--repair_subterminal_multirefs',
        action='store_true',
        required= False,
        help = 'Repair cases where a given path in the taxonomy contains sequences at subterminal nodes.'
    )

    parser.add_argument(
        '--repair_all_multirefs',
        action='store_true',
        required= False,
        help = 'Repair cases where a given path in the taxonomy contains sequences at any non-leaf nodes.'
    )

    return parser

def main():
    
    args = args_parser().parse_args()
    db_path = args.db_file or tempfile.NamedTemporaryFile(delete= not args.keep_db_file)
    
    kdb = KrakenDbBuilder(db_path= db_path)
    kdb.load_taxonomy_files(
        taxonomy_dir= args.taxonomy_path, 
        no_acc2taxid= args.no_acc2taxid
    )

    # Load the bulk of the reference genomes from FASTA files.
    # These are loaded into the "sequences" table of the sqlite DB without a category value, 
    for fasta_path in args.fasta_path:
        kdb.load_fasta_file(fasta_path)

    if args.filter_flu:
        kdb.filter_unnamed_unsegmented_flu()
        if args.filter_except:
            filter_except_list = args.filter_except
        else:
            filter_except_list = []
        kdb.filter_incomplete_flu(filter_except_patterns= filter_except_list)

    kdb.create_segmented_flu_taxonomy_nodes()
    kdb.assign_flu_taxonomy_nodes()

    multiref_paths, seen, multiref_data = kdb.find_multiref_paths(root_taxid = 10239)
    if not args.repair_subterminal_multirefs and not args.repair_all_multirefs:
        logging.info("Paths with sequences attached to non-leaf nodes will not be fixed.")
    if args.repair_subterminal_multirefs and not args.repair_all_multirefs:
        logging.info("Repairing only subterminal multiref cases.")
        kdb.repair_multiref_paths(multiref_paths, seen, "subterminal")
    if args.repair_all_multirefs:
        logging.info("Repairing all multiref cases.")
        kdb.repair_multiref_paths(multiref_paths, seen, "all")

    kdb.create_db_ready_dir(path = args.out_dir)

if __name__ == "__main__":
    exit(main())