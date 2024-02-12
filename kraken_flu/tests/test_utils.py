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
    
def test_FLU_REGEX():
    str = 'gi|59896308|gb|CY000029|Influenza A virus (A/New York/32/2003(H3N2)) segment 8, complete sequence'
    assert FLU_REGEX.search(str), 'regex matches expected pattern'
    
    str = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    assert FLU_REGEX.search(str), 'regex matches expected pattern'
    
def test_FLU_ISOLATE_NAME_REGEX():
    str = 'gi|59896308|gb|CY000029|Influenza A virus (A/New York/32/2003(H3N2)) segment 8, complete sequence'
    assert FLU_ISOLATE_NAME_REGEX.search(str), 'regex matches expected pattern'
    match = FLU_ISOLATE_NAME_REGEX.search(str)
    assert match.group(1) == 'A/New York/32/2003(H3N2)', 'regex extracts the correct pattern'
    
    str = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    assert FLU_ISOLATE_NAME_REGEX.search(str), 'regex matches expected pattern'
    match = FLU_ISOLATE_NAME_REGEX.search(str)
    assert match.group(1) == 'B/Texas/24/2020', 'regex extracts the correct pattern'
    
    str = '>gi|169731751|gb|CY030663|Influenza B virus (B/Tennessee/UR06-0431/2007) segment 1, complete sequence'
    assert FLU_ISOLATE_NAME_REGEX.search(str), 'regex matches expected pattern'
    match = FLU_ISOLATE_NAME_REGEX.search(str)
    assert match.group(1) == 'B/Tennessee/UR06-0431/2007', 'regex extracts the correct pattern'
    
    str = '>kraken:taxid|2697049|NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    assert FLU_ISOLATE_NAME_REGEX.search(str) is None, 'does not match'
    
def test_FLU_SEG_NUM_REGEX():
    str = 'gi|59896308|gb|CY000029|Influenza A virus (A/New York/32/2003(H3N2)) segment 8, complete sequence'
    assert FLU_SEG_NUM_REGEX.search(str), 'regex matches expected pattern'
    match = FLU_SEG_NUM_REGEX.search(str)
    assert match.group(1) == '8', 'regex extracts the correct pattern'
    
    str = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    assert FLU_SEG_NUM_REGEX.search(str), 'regex matches expected pattern'
    match = FLU_SEG_NUM_REGEX.search(str)
    assert match.group(1) == '2', 'regex extracts the correct pattern'
    
    str = '>kraken:taxid|2697049|NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    assert FLU_SEG_NUM_REGEX.search(str) is None, 'does not match'
    
    str = 'segment 1 of something that is not flu'
    assert FLU_SEG_NUM_REGEX.search(str) is None, 'does not match'
    
def test_KRAKEN_TAX_ID_REGEX():
    str = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    assert KRAKEN_TAX_ID_REGEX.search(str) is None, 'does not match'
    
    str = '>kraken:taxid|2697049|NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome'
    assert KRAKEN_TAX_ID_REGEX.search(str), 'regex matches expected pattern'
    match = KRAKEN_TAX_ID_REGEX.search(str)
    assert match.group(1) == '2697049', 'regex extracts the correct pattern'

def test_FLU_A_SUBTYPE_REGEX():
    assert FLU_A_SUBTYPE_PARTS_REGEX.search('H10N10')
    assert FLU_A_SUBTYPE_PARTS_REGEX.search('H18N10').group(1) == 'H18'
    assert FLU_A_SUBTYPE_PARTS_REGEX.search('H18N10').group(2) == 'N10'
    
def test_FLU_A_ISOLATE_NAME_REGEX():
    assert FLU_A_ISOLATE_NAME_REGEX.search('Influenza A virus (A/PR 8/34(H1N1))')
    match = FLU_A_ISOLATE_NAME_REGEX.search('Influenza A virus (A/PR 8/34(H1N1))')
    assert match.group(1) == 'A/PR 8/34(H1N1)'
    assert not FLU_A_ISOLATE_NAME_REGEX.search('Influenza B virus (STRAIN B/LEE/40)')
    
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
    
    match = FLU_DATA_REGEX.search('some segment virus')
    assert not match
    
def test_parse_flu():
    name = 'Influenza A virus (A/Puerto_Rico/8/34(H10N11))'
    flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( name )
    assert flu_type == 'A' ,'type'
    assert isolate_name == 'A/Puerto_Rico/8/34(H10N11)', 'isolate name'
    assert h_subtype == 'H10' , 'H subtype'
    assert n_subtype == 'N11' , 'N subtype'
    assert segment_number is None , 'segment number'
    
    name = 'Influenza B virus (STRAIN B/LEE/40) RNA 1'
    flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( name )
    assert flu_type == 'B' ,'type'
    assert isolate_name == 'STRAIN B/LEE/40', 'isolate name'
    assert h_subtype is None , 'H subtype'
    assert n_subtype is None , 'N subtype'
    assert segment_number == 1 , 'segment number'
    
    # should work with any input, including a full FASTA header like this
    name = '>gi|1834444346|gb|MT375832|Influenza B virus (B/Texas/24/2020) segment 2 polymerase PB2 (PB2) gene, complete cds'
    flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( name )
    assert flu_type == 'B' ,'type'
    assert isolate_name == 'B/Texas/24/2020', 'isolate name'
    assert h_subtype is None , 'H subtype'
    assert n_subtype is None , 'N subtype'
    assert segment_number == 2 , 'segment number'
    
    # THis is a rare case but it does happen: influenza without an isolate name
    # We do need to allow this so that we can perform an additional step later to try and fix this
    name = '>kraken:taxid|518987|NC_002204.1 Influenza B virus RNA 1, complete sequence'
    flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( name )
    assert flu_type == 'B' ,'type'
    assert isolate_name == None, 'no isolate name'
    assert h_subtype is None , 'H subtype'
    assert n_subtype is None , 'N subtype'
    assert segment_number == 1 , 'segment number'
    
