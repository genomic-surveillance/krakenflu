from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import subprocess
import logging 

from kraken_flu.src.utils import KRAKEN_TAX_ID_REGEX, NCBI_ACC_REGEX
from kraken_flu.src.utils import parse_flu
from kraken_flu.src.db import Db

"""
This module contains the functionality for loading sequence data from FASTA into the sqlite database    
"""

def load_fasta(db: Db, file_path:str, category:str=None, enforce_ncbi_acc:bool = False):
    """
    Main function. Orchestrates the loading of a FASTA sequence file into the database.

    Args:
        db: KrakenFlu::Db object, required
            This provides the connection to the DB
        
        file_path: str, required
            Path to the FASTA file we are loading into the DB
        
        category: str, optional
            If provided, will be used to set the column "category" in the "sequences" 
            table, which can be used later to create associations with taxonomy nodes.  
            This is used for cases where we load a specific FASTA file for a known virus 
            (type) and we want to save a hint in the DB for the taxonomy association later. 
            
        enforce_ncbi_acc: bool, optional, defaults to False
            If True, an exception is thrown if an NCBI acc ID cannot be found
            
    Returns:
        True on success
        
    Side effects:
        loads data into backend DB 
            
    """
    n_seqs = _get_num_records(file_path)
    logging.info( f'starting to upload {n_seqs} FASTA sequence records from {file_path} data to DB')
    if category:
        logging.info(f'setting category to "{category}"')
    
    with open( file_path, encoding="utf-8" ) as fh:
        with db.bulk_insert_buffer(table_name='sequences', buffer_size= 10000) as b:
            for record in SeqIO.parse(fh, "fasta"):
                header = record.description
                sequence = record.seq
                flu_type, ncbi_acc, kraken_taxid, is_flu, is_fluA, isolate_name, segment_number, h_subtype, n_subtype = _parse_header(header, enforce_ncbi_acc = enforce_ncbi_acc)
                
                # For the DB, we just want to store the integer of the H and N subtype
                h_subtype = h_subtype and h_subtype.replace('H', '')
                n_subtype = n_subtype and n_subtype.replace('N', '')
                sequence = str(sequence)
                seq_len = len(sequence)
                percent_n_bases = _calculate_percent_n( sequence= sequence)
                
                n_inserted = b.add_row(
                    {
                        'fasta_header': header,
                        'dna_sequence': sequence,
                        'seq_length': seq_len,
                        'percent_n': percent_n_bases,
                        'category': category,
                        'flu_type': flu_type,
                        'ncbi_acc': ncbi_acc,
                        'original_tax_id': kraken_taxid,
                        'is_flu': int(is_flu),
                        'flu_name': isolate_name,
                        'segment_number': segment_number, 
                        'flu_a_h_subtype': h_subtype,
                        'flu_a_n_subtype': n_subtype,
                        'include':int(True)
                    }
                )
                if n_inserted > 0:
                    logging.info(f'flushed {n_inserted} records to DB')

    logging.info( f'finished uploading sequence records to DB')
    return True

def _calculate_percent_n(sequence:str):
    """
    Calculates the percentage of N bases in a sequence.  

    Args:
        sequence: str, required
            The sequence string
    """
    sequence_length = len(sequence)
    n_count = sequence.upper().count('N')
    return (n_count / sequence_length) * 100

def _parse_header(header:str, enforce_ncbi_acc:bool= False):
    """
    Parse a FASTA header
    
    Parameters:
        header: str, required
            the FASTA header string
            
        enforce_ncbi_acc: bool, optional, defaults to False
            If True, an exception is thrown if an NCBI acc ID cannot be found
            
    Returns tuple of:
        ncbi_acc (str), 
        kraken_taxid (None or int), 
        is_flu (bool),
        is_fluA (bool),
        flu_isolate_name (str), 
        flu_segment_number (None or int),
        fluA_H_subtype
        fluA_N_subtype
            
    """
    try:
        ncbi_acc = NCBI_ACC_REGEX.search( header ).group(0)
        ncbi_acc = re.sub(r'^gb\|', '', ncbi_acc)
    except AttributeError:
        if enforce_ncbi_acc:
            raise ValueError(f"could not parse NCBI acc ID from FASTA header '{header}' - use enforce_ncbi_acc=False to ignore this error")
        else:
            ncbi_acc = None

    if (match := KRAKEN_TAX_ID_REGEX.search( header )) is not None:
        kraken_taxid = int(match.group(1))
    else:
        kraken_taxid = None
        
    is_flu, flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( header )
    if flu_type and flu_type == 'A':
            is_fluA = True
    else:
        is_fluA = False
        
    return flu_type, ncbi_acc, kraken_taxid, is_flu, is_fluA, isolate_name, segment_number, h_subtype, n_subtype

def _get_num_records(file_path):
    """
    Counts FASTA records in file file_path
    """
    p = subprocess.run(f"grep -c '^>' {file_path}", shell=True, check=True, capture_output=True, encoding='utf-8')
    if not p.returncode:
        return int(p.stdout)
    else:
        raise Exception("failed to runcommand to count FASTA records in file")