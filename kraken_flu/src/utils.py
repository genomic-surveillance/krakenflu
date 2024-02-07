import re

# regular expressions for FASTA header and taxonomy names parsing
FLU_REGEX = re.compile(r'Influenza[ _][AB].+')
FLU_ISOLATE_NAME_REGEX = re.compile(r'Influenza[ _][AB].*?\(([A-Za-z\-_ /[0-9]*?(\(H[0-9]+N[0-9]+\))?)\)')
FLU_A_ISOLATE_NAME_REGEX = re.compile(r'Influenza[ _]A.*?\((A\/[A-Za-z\-_ /[0-9]*?(\(H[0-9]+N[0-9]+\))?)\)')
FLU_SEG_NUM_REGEX = re.compile(r'Influenza[ _][AB].+ (?:segment|RNA) ([1-8])')
KRAKEN_TAX_ID_REGEX = re.compile(r'kraken:taxid\|([0-9]+)\|')
# use GenBank (gb) number or RefSeq accession such as NC_xxxxxx (RefSeq chromosome)
# could potentially be stricter with the RefSeq IDs as we are only interested in NC_xxxxx(?)
NCBI_ACC_REGEX = re.compile(r'gb\|[A-Z]+_?[0-9]+|[A-Z]{2}_[0-9]{6,}\.[0-9]')
FLU_A_SUBTYPE_PARTS_REGEX = re.compile(r'(H[0-9]+)(N[0-9]+)')

# TODO: this Regex should eventually replace all of the above FLU regexes (combines all in one)
FLU_DATA_REGEX = re.compile(
    # example: Influenza A virus (A/Puerto Rico/8/1934(H1N1)) segment 1
    r'Influenza[ _]([A-D])'+ # influenza type, capture group 1: A
    r'.*?\('+ # most likely "virus ("
        r'([A-Za-z\-_ /[0-9]*?'+ # isolate name, capture group 2: A/Puerto Rico/8/1934(H1N1)
            r'(?:'+ # optional non capturing group for HxNx subtype (only present in flu A)
                r'\((H[0-9]+)'+ # H subtype, capture group 3: H1
                r'(N[0-9]+)\)'+ # N subtype, capture group 4: N1
            r')?'+ # end of optional non-capturing group for subtype
        r')'+ # end of capture group 2 (isolate name)
    r'\)'+ # isolate name closing parenthesis
    r'(?:.*(?:segment|RNA) ([1-8]))?' ) # segment number, capture group 5 (optional)

def parse_flu( in_str:str ):
    """
    Parses key information about flu isolates from a string, which could be a FASTA
    header or a name in the taxonomy.
    
    Parameters:
        in_str: str, required
            The string to parse
            
    Returns:
        flu_type, isolate name, H subtype, N subtype, segment number
        
        where:
            flu_type: str ('A', 'B', ...)
            H subtype: str ('H1', 'H2', 'H3', ...)
            N subtype: str ('N1', 'N2', ...)
            isolate name: str ('A/PR8/1934(H1N1)', ... )
            segment number: int
            
            Each return value can be None if the pattern does not match
    """
    flu_type = h_subtype = n_subtype = isolate_name = segment_number = None
    match = FLU_DATA_REGEX.match( in_str )
    if match:
        flu_type = match.group( 1 )
        isolate_name = match.group( 2 )
        h_subtype = match.group( 3 )
        n_subtype = match.group( 4 )
        segment_number = match.group( 5 )
        if segment_number:
            segment_number = int(segment_number)
    return flu_type, isolate_name, h_subtype, n_subtype, segment_number