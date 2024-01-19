from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import re
import argparse
import random

"""
From a fresh Kraken2 DB preparation directory, extract the H1N1, H3N2, RSV and COVID
sequences into a new file from the existing library.fa.
This script is here because it is ONLY used to create the test fixtures file.
It was run like this to create the test file:
make_test_fasta.py --in_fasta /lustre/scratch126/gsu/team112/personal/fs5/rvi_dev/krakenDBs/customDB-viralrefseq_replaced_flu_segments/library/viral/library.fna --out_fasta library/viral/library.fna

followed by:
grep ">" library/viral/library.fna | cut -d" " -f 1 | cut -d"|" -f 3 > ~/tmp/nc_ids
grep -f ~/tmp/nc_ids /lustre/scratch126/gsu/team112/personal/fs5/rvi_dev/krakenDBs/customDB-viralrefseq_replaced_flu_segments/library/viral/prelim_map.txt > library/viral/prelim_map.txt

make a subset of taxonomy files names.dmp and nodes.dmp
the regex allows some extra records to 'slip through' but that doesn't matter for these files 
cut -d"|" -f 2 library/viral/prelim_map.txt | sort | uniq > ~/tmp/tax_ids
grep -w -f ~/tmp/tax_ids /lustre/scratch126/gsu/team112/personal/fs5/rvi_dev/krakenDBs/customDB-viralrefseq_replaced_flu_segments/taxonomy/names.dmp > taxonomy/names.dmp
grep -w -f ~/tmp/tax_ids /lustre/scratch126/gsu/team112/personal/fs5/rvi_dev/krakenDBs/customDB-viralrefseq_replaced_flu_segments/taxonomy/nodes.dmp > taxonomy/nodes.dmp

add the parent taxa
grep -e 'H1N1 subtype' /lustre/scratch126/gsu/team112/personal/fs5/rvi_dev/krakenDBs/customDB-viralrefseq_replaced_flu_segments/taxonomy/names.dmp >> tests/fixtures/kraken_ncbi_data/taxonomy/names.dmp 
(devenv) fs5@farm5-head1:kraken_db_maker$ grep -e 'H3N2 subtype' /lustre/scratch126/gsu/team112/personal/fs5/rvi_dev/krakenDBs/customDB-vir
alrefseq_replaced_flu_segments/taxonomy/names.dmp >> tests/fixtures/kraken_ncbi_data/taxonomy/names.dmp 
(devenv) fs5@farm5-head1:kraken_db_maker$ grep -P '\tInfluenza B virus\t' /lustre/scratch126/gsu/team112/personal/fs5/rvi_dev/krakenDBs/customDB-viralrefseq_replaced_flu_segments/taxonomy/names.dmp >> tests/fixtures/kraken_ncbi_data/taxonomy/names.dmp


It is also possible to downsample the input fasta.
This command will keep ~10% of the sequences

make_test_fasta.py --in_fasta in.fasrta --out_fasta out.fasta --downsample 0.1
"""

REGEXES = [
    re.compile(r'Influenza A.*H1N1'),
    re.compile(r'Influenza A.*H3N2'),
    re.compile(r'Influenza B'),
    re.compile(r'NC_001803\.1'),
    re.compile(r'NC_001781\.1'),
    re.compile(r'NC_045512\.2') ]

def args_parser():
        
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '--in_fasta',
        action = 'store',
        required = True,
        metavar = 'FILE', 
        type = str,
        help='multi-sequence FASTA to process')
    
    parser.add_argument(
        '--out_fasta',
        action = 'store',
        required = True,
        metavar = 'FILE', 
        type = str,
        help='output file for matching sequences')
    
    parser.add_argument(
        '--downsample',
        action = 'store',
        required = False,
        type = float,
        help = 'proportion (0-1) of the sequences to keep. Sequences will be chosen randomly'
    )
    
    return parser

def main():

    args = args_parser().parse_args()
    
    with open( args.in_fasta ) as in_handle:
        with open( args.out_fasta, 'w' ) as out_handle:
            for record in SeqIO.parse(in_handle, "fasta"):
                for regex in REGEXES:
                    if regex.search( record.description):
                        if not args.downsample or random.random() <= args.downsample: 
                            SeqIO.write(record, out_handle, "fasta")
                      
if __name__ == "__main__":
    exit(main())
