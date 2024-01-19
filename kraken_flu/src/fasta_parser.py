from Bio import SeqIO
import re
import os.path
from cached_property import cached_property
import logging 

logging.basicConfig( format='%(asctime)s %(message)s', level=logging.DEBUG )

class FastaParser():   
    """
    This class parsers FASTA and extracts the following information from the headers:
        - pre-assigned kraken2 taxid, if present
        - NCBI accession (Genbank unique ID)
        - boolean: is_flu: True if the header indicates an influenza sequence
        - flu_name: populated in case it is an influenza name
            examples: A/New York/392/2004(H3N2), B/Houston/B850/2005
        - flu_segment: populated for influenza only
        - length of the sequence (the sequence itself is not stored)
    
    Parameters:
        fasta_file_path: str, required
            path to the FASTA file that is to be filtered
        
    """
    
    def __init__( self, fasta_file_path: str ):
        if not os.path.isfile( fasta_file_path ):
            raise ValueError(f'file { fasta_file_path } does not exist or is not a file')
        self.fasta_file_path = fasta_file_path
