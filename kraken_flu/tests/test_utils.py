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
