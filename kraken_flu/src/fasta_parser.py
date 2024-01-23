from Bio import SeqIO
import re
import os.path
from cached_property import cached_property
import logging 

logging.basicConfig( format='%(asctime)s %(message)s', level=logging.DEBUG )

class FastaParser():   
    """
    This class parses FASTA and extracts the following information from the headers:
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
        
    # regular expressions for FASTA header parsing
    FLU_REGEX = re.compile(r'Influenza[ _][AB].+\(.+\)')
    FLU_ISOLATE_NAME_REGEX = re.compile(r'Influenza[ _][AB].*?\(([A-Za-z\-_ /[0-9]*?(\(H[0-9]+N[0-9]+\))?)\)')
    FLU_SEG_NUM_REGEX = re.compile(r'Influenza[ _][AB].+ segment ([1-8])')
    KRAKEN_TAX_ID_REGEX = re.compile(r'kraken:taxid\|([0-9]+)\|')
    # use GenBank (gb) number or RefSeq accession such as NC_xxxxxx (RefSeq chromosome)
    # could potentially be stricter with the RefSeq IDs as we are only interested in NC_xxxxx(?)
    NCBI_ACC_REGEX = re.compile(r'gb\|[A-Z]+_?[0-9]+|[A-Z]{2}_[0-9]{6,}\.[0-9]')


    def __init__( self, fasta_file_path: str ):
        if not os.path.isfile( fasta_file_path ):
            raise ValueError(f'file { fasta_file_path } does not exist or is not a file')
        self.fasta_file_path = fasta_file_path

    def _parse_header( self, header):
        """
        Parse a FASTA header
        
        Parameters:
            header: str, required
                the FASTA header string
                
        Returns tuple of:
            ncbi_acc (str), kraken_taxid (None or int), is_flu (bool), flu_isolate_name (str), flu_segment_number (None or int) 
                
        """
        try:
            ncbi_acc = self.NCBI_ACC_REGEX.search( header ).group(0)
            ncbi_acc = re.sub(r'^gb\|', '', ncbi_acc)
        except AttributeError:
            # an NCBI accession ID is required
            raise ValueError(f"could not parse NCBI acc ID from FASTA header '{header}'")

        if (match := self.KRAKEN_TAX_ID_REGEX.search( header )) is not None:
            kraken_taxid = int(match.group(1))
        else:
            kraken_taxid = None
            
        flu_isolate_name = None
        flu_segment_number = None
        if self.FLU_REGEX.search( header ):
            is_flu = True
            if (match := self.FLU_SEG_NUM_REGEX.search( header )) is not None:
                flu_segment_number = int(match.group(1))
            if (match := self.FLU_ISOLATE_NAME_REGEX.search( header )) is not None:
                flu_isolate_name = match.group(1)
        else:
            is_flu = False
            
        return ncbi_acc, kraken_taxid, is_flu, flu_isolate_name, flu_segment_number 

    @cached_property
    def header_data( self ):
        """
        Returns FASTA data as a list of dicts, each dict has the following keys:
            - orig_head: the original un-modified FASTA header
            - mod_head: modified header:
                - any kraken tax ID removed
                - NCBI accession ID preserved
                - flu genomes renamed to: ISOLATE_NAME SEGMENT_NUMBER 
            - ncbi_acc: NCBI accession ID
            - taxid: kraken taxid if present
            - seqlen: length of the sequence
            The following are only applied to flu sequences (None otherwise):
            - is_flu: True if this is a flu sequence
            - flu_name: flu isolate name such as 'A/New York/32/2003(H3N2)' or 'B/Texas/24/2020'
            - flu_seg_num: segment number 
        
        Returns:
            list of dicts, see above
        
        """
        logging.info( f'parsing FASTA headers from { self.fasta_file_path }')
        data = []
        n_all = 0
        n_flu = 0
        with open( self.fasta_file_path ) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                orig_header = record.description
                seqlen = len(record.seq)
                ncbi_acc, kraken_taxid, is_flu, flu_isolate_name, flu_segment_number = self._parse_header(orig_header)
                
                n_all+=1
                if is_flu:
                    n_flu+=1
                    mod_header = ' '.join( 
                        ['gb|'+ ncbi_acc +'|', 
                        'Influenza', 
                        flu_isolate_name, 
                        'segment', 
                        str(flu_segment_number)])
                else:
                    mod_header = self.KRAKEN_TAX_ID_REGEX.sub('', orig_header)

                data.append({
                    'orig_head': orig_header,
                    'mod_head': mod_header,
                    'seqlen': seqlen,
                    'ncbi_acc': ncbi_acc,
                    'is_flu': is_flu,
                    'taxid': kraken_taxid,
                    'flu_name': flu_isolate_name,
                    'flu_seg_num': flu_segment_number })
                
                
        
        logging.info( f'found {n_all} sequences, {n_flu} of which are influenza')
        return data
    
    