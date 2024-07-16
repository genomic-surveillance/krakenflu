from kraken_flu.src.utils import *

def test_NCBI_ACC_REGEX():
    str = 'gi|59896308|gb|CY000029|Influenza A virus (A/New York/32/2003(H3N2)) segment 8, complete sequence'
    assert NCBI_ACC_REGEX.search(str), 'regex matches expected pattern'
    match = NCBI_ACC_REGEX.search(str)
    assert match.group(0) == 'gb|CY000029', 'regex extracts the correct pattern'
    
    str = '>kraken:taxid|2697049|NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    assert NCBI_ACC_REGEX.search(str), 'regex matches expected pattern'
    match = NCBI_ACC_REGEX.search(str)
    assert match.group(0) == 'NC_045512.2', 'regex extracts the correct pattern'
    
def test_KRAKEN_TAX_ID_REGEX():
    str = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    assert KRAKEN_TAX_ID_REGEX.search(str) is None, 'does not match'
    
    str = '>kraken:taxid|2697049|NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    assert KRAKEN_TAX_ID_REGEX.search(str), 'regex matches expected pattern'
    match = KRAKEN_TAX_ID_REGEX.search(str)
    assert match.group(1) == '2697049', 'regex extracts the correct pattern'
    
def test_FLU_DATA_REGEX():
    match = FLU_DATA_REGEX.search('Influenza A virus (A/PR 8/34(H1N1)) segment 2')
    assert match
    assert match.group(1) == 'A', 'type'
    assert match.group(2) == 'A/PR 8/34(H1N1)', 'isolate name'
    assert match.group(3) == 'H1', 'H subtype'
    assert match.group(4) == 'N1', 'N subtype'
    assert match.group(5) == '2', 'segment number'
    
    match = FLU_DATA_REGEX.search('Influenza A virus (A/Puerto_Rico/8/34(H10N11))')
    assert match
    assert match.group(1) == 'A', 'type'
    assert match.group(2) == 'A/Puerto_Rico/8/34(H10N11)', 'isolate name'
    assert match.group(3) == 'H10', 'H subtype'
    assert match.group(4) == 'N11', 'N subtype'
    assert match.group(5) == None, 'segment number'
    
    match = FLU_DATA_REGEX.search('Influenza B virus (STRAIN B/LEE/40) RNA 1')
    assert match
    assert match.group(1) == 'B', 'type'
    assert match.group(2) == 'STRAIN B/LEE/40', 'isolate name'
    assert match.group(3) == None, 'H subtype'
    assert match.group(4) == None, 'N subtype'
    assert match.group(5) == '1', 'segment number'

    # this caused an issue because literal dot wasn't allowed in FLU_DATA_REGEX        
    match = FLU_DATA_REGEX.search('Influenza A virus (A/swine/Thailand/KU5.1/2004(H3N2))')
    assert match
    assert match.group(1) == 'A', 'type'
    assert match.group(2) == 'A/swine/Thailand/KU5.1/2004(H3N2)', 'isolate name'
    assert match.group(3) == 'H3', 'H subtype'
    assert match.group(4) == 'N2', 'N subtype'
    assert match.group(5) == None, 'segment number'
    
    match = FLU_DATA_REGEX.search('Influenza A virus (A/Fiji/15899/83)')
    assert match
    assert match.group(1) == 'A', 'type'
    assert match.group(2) == 'A/Fiji/15899/83', 'isolate name'
    assert match.group(3) == None, 'H subtype'
    assert match.group(4) == None, 'N subtype'
    assert match.group(5) == None, 'segment number'
    
    match = FLU_DATA_REGEX.search('Influenza A virus (A/Hawaii/21/2001(N2))')
    assert match
    assert match.group(1) == 'A', 'type'
    assert match.group(2) == 'A/Hawaii/21/2001(N2)', 'isolate name'
    assert match.group(3) == None, 'H subtype'
    assert match.group(4) == 'N2', 'N subtype'
    assert match.group(5) == None, 'segment number'    
    
    match = FLU_DATA_REGEX.search('Influenza A virus (A/Iowa/10/2015(H3))')
    assert match
    assert match.group(1) == 'A', 'type'
    assert match.group(2) == 'A/Iowa/10/2015(H3)', 'isolate name'
    assert match.group(3) == 'H3', 'H subtype'
    assert match.group(4) == None, 'N subtype'
    assert match.group(5) == None, 'segment number'

    match = FLU_DATA_REGEX.search('Influenza A virus (A/swine/Argentina/CIP112-C93.85/2014(H1N1))')
    assert match
    assert match.group(1) == 'A', 'type'
    assert match.group(2) == 'A/swine/Argentina/CIP112-C93.85/2014(H1N1)', 'isolate name'
    assert match.group(3) == 'H1', 'H subtype'
    assert match.group(4) == 'N1', 'N subtype'
    assert match.group(5) == None, 'segment number'
    
    match = FLU_DATA_REGEX.search('Influenza A virus (A/swine/Iowa/A01477719/2014(mixed))')
    assert match
    assert match.group(1) == 'A', 'type'
    assert match.group(2) == 'A/swine/Iowa/A01477719/2014(mixed)', 'isolate name'
    assert match.group(3) == None, 'H subtype'
    assert match.group(4) == None, 'N subtype'
    assert match.group(5) == None, 'segment number'

    match = FLU_DATA_REGEX.search('some segment virus')
    assert not match

    
def test_parse_flu():
    name = 'Influenza A virus (A/Puerto_Rico/8/34(H10N11))'
    is_flu, flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( name )
    assert flu_type == 'A' ,'type'
    assert isolate_name == 'A/Puerto_Rico/8/34(H10N11)', 'isolate name'
    assert h_subtype == 'H10' , 'H subtype'
    assert n_subtype == 'N11' , 'N subtype'
    assert segment_number is None , 'segment number'
    
    name = 'Influenza B virus (STRAIN B/LEE/40) RNA 1'
    is_flu, flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( name )
    assert flu_type == 'B' ,'type'
    assert isolate_name == 'STRAIN B/LEE/40', 'isolate name'
    assert h_subtype is None , 'H subtype'
    assert n_subtype is None , 'N subtype'
    assert segment_number == 1 , 'segment number'
    
    # should work with any input, including a full FASTA header like this
    name = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    is_flu, flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( name )
    assert flu_type == 'B' ,'type'
    assert isolate_name == 'B/Texas/24/2020', 'isolate name'
    assert h_subtype is None , 'H subtype'
    assert n_subtype is None , 'N subtype'
    assert segment_number == 2 , 'segment number'
    
    # THis is a rare case but it does happen: influenza without an isolate name
    # We do need to allow this so that we can perform an additional step later to try and fix this
    name = '>kraken:taxid|518987|NC_002204.1 Influenza B virus RNA 1, complete sequence'
    is_flu, flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( name )
    assert flu_type == 'B' ,'type'
    assert isolate_name == None, 'no isolate name'
    assert h_subtype is None , 'H subtype'
    assert n_subtype is None , 'N subtype'
    assert segment_number == 1 , 'segment number'
    
    # Flu B without a segment number
    name = 'Influenza B virus (B/Lee/1940)'
    flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( name )
    assert flu_type == 'B' ,'type'
    assert isolate_name == 'B/Lee/1940', 'isolate name'
    assert h_subtype is None , 'H subtype'
    assert n_subtype is None , 'N subtype'
    assert segment_number == None , 'segment number'
    
    # Real-world case of a genome segment that does not have a segment number but 
    # does have a gene name that can be translated into a number
    name = '>NC_007360.1 Influenza A virus (A/Goose/Guangdong/1/96(H5N1)) nucleocapsid protein (NP) gene, complete cds'
    is_flu, flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( name )
    assert flu_type == 'A' ,'type'
    assert isolate_name == 'A/Goose/Guangdong/1/96(H5N1)', 'isolate name'
    assert h_subtype == 'H5' , 'H subtype'
    assert n_subtype == 'N1' , 'N subtype'
    assert segment_number == 5 , 'segment number correctly derived from gene name'
    
    # additional tests for isolate names that caused issues
    name = 'Influenza A virus (A/Fiji/15899/83)'
    is_flu, flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( name )
    assert flu_type == 'A' ,'type'
    assert isolate_name == 'A/Fiji/15899/83', 'isolate name'
    assert h_subtype == None , 'H subtype'
    assert n_subtype == None , 'N subtype'
    assert segment_number == None , 'segment number'
    
    name = 'Influenza A virus (A/Hawaii/21/2001(N2))'
    is_flu, flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( name )
    assert flu_type == 'A' ,'type'
    assert isolate_name == 'A/Hawaii/21/2001(N2)', 'isolate name'
    assert h_subtype == None , 'H subtype'
    assert n_subtype == 'N2' , 'N subtype'
    assert segment_number == None , 'segment number'
    
    name = 'Influenza A virus (A/swine/Thailand/KU5.1/2004(H3N2))'
    is_flu, flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( name )
    assert flu_type == 'A' ,'type'
    assert isolate_name == 'A/swine/Thailand/KU5.1/2004(H3N2)', 'isolate name'
    assert h_subtype == 'H3' , 'H subtype'
    assert n_subtype == 'N2' , 'N subtype'
    assert segment_number == None , 'segment number'