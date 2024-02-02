import re

# regular expressions for FASTA header and taxonomy names parsing
FLU_REGEX = re.compile(r'Influenza[ _][AB].+')
FLU_ISOLATE_NAME_REGEX = re.compile(r'Influenza[ _][AB].*?\(([A-Za-z\-_ /[0-9]*?(\(H[0-9]+N[0-9]+\))?)\)')
FLU_SEG_NUM_REGEX = re.compile(r'Influenza[ _][AB].+ (?:segment|RNA) ([1-8])')
KRAKEN_TAX_ID_REGEX = re.compile(r'kraken:taxid\|([0-9]+)\|')
# use GenBank (gb) number or RefSeq accession such as NC_xxxxxx (RefSeq chromosome)
# could potentially be stricter with the RefSeq IDs as we are only interested in NC_xxxxx(?)
NCBI_ACC_REGEX = re.compile(r'gb\|[A-Z]+_?[0-9]+|[A-Z]{2}_[0-9]{6,}\.[0-9]')
FLU_A_SUBTYPE_REGEX = re.compile(r'H([0-9]+)N([0-9]+)')