from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import os.path
import logging 

from kraken_flu.src.utils import KRAKEN_TAX_ID_REGEX, NCBI_ACC_REGEX
from kraken_flu.src.utils import parse_flu
from kraken_flu.src.db import Db

def load_fasta(db: Db, file_path:str, category:str=None):
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
            
    Returns:
        True on success
        
    Side effects:
        loads data into backend DB 
            
    """
    
    logging.info( f'uploading FASTA to DB from { file_path }')
    n=0
    with open( file_path ) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            n+=1
            header = record.description
            sequence = record.seq
            seqlen = len( sequence )
            flu_type, ncbi_acc, kraken_taxid, is_flu, is_fluA, isolate_name, segment_number, h_subtype, n_subtype = _parse_header(header)
            
            # For the DB, we jsut want to store the integer of the H and N subtype
            h_subtype = h_subtype and h_subtype.replace('H', '')
            n_subtype = n_subtype and n_subtype.replace('N', '')
            
            db.add_sequence(
                fasta_header= header,
                dna_sequence= str(sequence),
                category= category,
                flu_type= flu_type,
                ncbi_acc= ncbi_acc,
                original_taxid= kraken_taxid,
                is_flu= is_flu,
                isolate_name= isolate_name,
                segment_number= segment_number, 
                h_subtype= h_subtype,
                n_subtype= n_subtype
            )

    logging.info( f'finished uploading {n} sequences from { file_path }')
    
    return True

def _parse_header(header:str):
    """
    Parse a FASTA header
    
    Parameters:
        header: str, required
            the FASTA header string
            
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
        # an NCBI accession ID is required
        raise ValueError(f"could not parse NCBI acc ID from FASTA header '{header}'")

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

