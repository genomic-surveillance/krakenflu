from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re

from kraken_flu.src.db import Db

def write_fasta( db:Db, path:str, included_only:bool=True ):
        """
        Writes database content from sequences table to a new FASTA file.
        Only records with the "include" flag will be included. 
        
        FASTA header: if a modified FASTA header was recorded (which we do for influenza), the 
        modified version is used. Otherwise, the original header is used except that any existing 
        kraken:taxid tag is removed and replaced with the tax_id assigned to the sequence, if one exists.  
        Where no tax_id was assigned to the sequence but the input FASTA did have a kraken:taxid tag, the  
        latter is used. 
        
        Paremeters:
            db: Db, required
                KrakenFlu::Db object
                
            path: str, required
                Path to the file we are creating
                
            included_only: bool, optional, defaults to True
                If True, ony sequences with the "include=1" flag will be writte to
                the output file
                
        Returns:
            True on success
            
        Side-effects:
            Writes to file
            
        """
        seq_it = db.all_sequences_iterator(included_only= included_only)
        with open( path, 'w' ) as out_fh:
            for seq in seq_it:
                # use modified header if we have one
                fasta_header = seq['mod_fasta_header'] or seq['fasta_header']
                fasta_header = _remove_taxid(fasta_header)
                taxid = seq['tax_id'] or seq['original_tax_id']
                if taxid:
                    fasta_header = _add_taxid(fasta_header, taxid)
                sr = SeqRecord(
                    Seq(seq['dna_sequence']),
                    id= fasta_header,
                    description= ''
                )
                SeqIO.write( sr, out_fh, "fasta" )

        return True
    
def _remove_taxid(header):
    """
    Remove any existing kraken:taxid tag from a fasta header
    """
    return re.sub(r'kraken:taxid\|[0-9]+\|','', header)

def _add_taxid(header, taxid):
    """
    Add a kraken:taxid tag to the fasta header
    """
    return 'kraken:taxid|' + str(taxid) + '|' + header