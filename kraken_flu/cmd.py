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
        '--filter_flu_taxonomy',
        action = 'store_true',
        help = 'apply the influenza taxonomy filter: remove sequences that bypass the custom flu taxonomy. Must be used with do_complete_linkage'        
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
        '--do_full_linkage',
        action = 'store_true',
        help = 'run linkage of all sequences via NCBI accession ID to NCBI taxon ID for sequences that are not linked otherwise'        
    )

    parser.add_argument(
        '--prune_db',
        action = 'store_true',
        help = 'remove unnecessary data from the kraken-flu backend DB at the end of the process'        
    )

    # RSV options:
    # These three options depend on each other and have to be used in conjunction but there
    # is no argparse functionality for this as of now, so the logic is implemented in "main"
    parser.add_argument(
        '--rsv_a_sequences',
        action = 'store',
        required = False,
        metavar= 'FILE',
        type= str,
        help = 'file of known RSV A sequences. Must be used together with --rsv_b_sequences'
    )
    
    parser.add_argument(
        '--rsv_b_sequences',
        action = 'store',
        required = False,
        metavar= 'FILE',
        type= str,
        help = 'file of known RSV B sequences. Must be used together with --rsv_a_sequences'
    )

    parser.add_argument(
        '--rsv_size_filter',
        action = 'store_true',
        help = 'Must be used together with --rsv_a/b_sequences. Filters sequences in files for genome completeness (size)'        
    )
    
    parser.add_argument(
        '--max_percent_n',
        action = 'store',
        required = False,
        type= float,
        help = 'cutoff (in percent of bases) for maximum allowed N bases in any sequence'
    )
    
    parser.add_argument(
        '--deduplicate',
        action = 'store_true',
        help = 'de-duplicate sequences and keep only a single record for every group'        
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

    parser.add_argument(
        '--multiref_root_tax_id',
        type = int,
        default = 10239, #viruses
        required= False,
        help = 'The taxid expected to be in all paths of interest. Only paths containing this tax_id will be examined for multi-references and fixed. Defaults to 10239 (the taxid for "Viruses")'
    )

    return parser

def main():
    
    args = args_parser().parse_args()
    db_path = args.db_file or tempfile.NamedTemporaryFile(delete= not args.keep_db_file)
    
    # sanity checks for arguments
    if args.rsv_size_filter:
        if not args.rsv_a_sequences or not args.rsv_a_sequences:
            raise ValueError('Option --rsv_size_filter must be used together with --rsv_a_sequence and --rsv_b_sequence')
    if (args.rsv_a_sequences and not args.rsv_b_sequences) or (args.rsv_b_sequences and not args.rsv_a_sequences):
        raise ValueError('parameters --rsv_a_sequences and --rsv_b_sequences must be used together')
    
    if args.filter_flu_taxonomy and not args.do_full_linkage:
        raise ValueError('option --filter_flu_taxonomy can only be used with --do_full_linkage')
    
    # TODO: move the rest into a single function in KrakenDbBuilder, perhaps called "build"
    # We need to run some integration tests and this would be difficult in the current setup
    
    kdb = KrakenDbBuilder(db_path= db_path)
    kdb.load_taxonomy_files(
        taxonomy_dir= args.taxonomy_path, 
        no_acc2taxid= args.no_acc2taxid
    )

    # Load the bulk of the reference genomes from FASTA files.
    # These are loaded into the "sequences" table of the sqlite DB without a category value, 
    for fasta_path in args.fasta_path:
        kdb.load_fasta_file(file_path= fasta_path,  enforce_ncbi_acc= True)
        
    if args.filter_flu:
        kdb.filter_unnamed_unsegmented_flu()
        if args.filter_except:
            filter_except_list = args.filter_except
        else:
            filter_except_list = []
        kdb.filter_incomplete_flu(filter_except_patterns= filter_except_list)
        kdb.filter_flu_a_wo_subtype()
        
    kdb.create_segmented_flu_taxonomy_nodes()
    kdb.assign_flu_taxonomy_nodes()

    if args.rsv_a_sequences and args.rsv_b_sequences:
        kdb.load_fasta_file(file_path= args.rsv_a_sequences, category= 'RSV A', enforce_ncbi_acc= False)
        kdb.load_fasta_file(file_path= args.rsv_b_sequences, category= 'RSV B', enforce_ncbi_acc= False)
        kdb.create_rsv_taxonomy( rsv_size_filter= args.rsv_size_filter )
    
    if args.deduplicate:
        kdb.deduplicate_sequences()
        
    if args.max_percent_n:
        kdb.filter_max_percent_n(args.max_percent_n)
    
    if args.do_full_linkage:
        kdb.link_all_unlinked_sequences_to_taxonomy_nodes()
        
        # If we are creating a custom RSV tree and doing the full linkage, we can end up with 
        # RefSeq sequences being linked again to high-level RSV taxonomy nodes via the Genabk acc lookup
        # and the taxonomy assigned to the sequence on NCBI.  
        # THis call will remove all RSV sequences linked to higher-order RSV taxa except for those linked 
        # to hRSV A/B
        if args.rsv_a_sequences or args.rsv_b_sequences:
            kdb.filter_out_sequences_linked_to_high_level_rsv_nodes()
            
        # as for RSV: filter out flu sequences that are linked to high-level taxa (not our custom 
        # taxonomy nodes). These will include sequences that cannot even be recognised as flu from the 
        # name (don't contain keywords)
        if args.filter_flu_taxonomy:
            kdb.filter_out_sequences_linked_to_high_level_flu_nodes()
    
    multiref_paths, seen, multiref_data = kdb.find_multiref_paths(root_taxid = args.multiref_root_tax_id)
    if not args.repair_subterminal_multirefs and not args.repair_all_multirefs:
        logging.info("Paths with sequences attached to non-leaf nodes will not be fixed.")
    if args.repair_subterminal_multirefs and not args.repair_all_multirefs:
        logging.info("Repairing only subterminal multiref cases.")
        kdb.repair_multiref_paths(multiref_paths, seen, "subterminal")
    if args.repair_all_multirefs:
        logging.info("Repairing all multiref cases.")
        kdb.repair_multiref_paths(multiref_paths, seen, "all")

    kdb.create_db_ready_dir(path = args.out_dir)
    
    # if a DB prune is requested, carry it out but only if we are actually keeping the DB
    # because it makes no sense to do this sort of housekeeping on a DB that is deleted anyway
    if args.prune_db and args.keep_db_file:
        kdb.prune_db()

if __name__ == "__main__":
    exit(main())