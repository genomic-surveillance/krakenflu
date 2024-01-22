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
        
    { genome_name: 
        [
            'NCBI_acc': NCBI accession ID,  
        {
            
            
            'kraken_taxid': taxon ID if present, otherwise None,
        }
    
    Parameters:
        fasta_file_path: str, required
            path to the FASTA file that is to be filtered
        
    """
        
    # These regexes define which FASTA header patterns we treat as flu genomes that need to 
    # be split into one taxon per segment
    # Flu genomes that are not covered by this pattern will still be present in the final
    # genomes but they will not be modified, so will not be split into segments in the kraken2 report
    FLU_REGEX = re.compile(r'Influenza [AB].+\(.+\)')
    FLU_ISOLATE_NAME_REGEX = re.compile(r'Influenza [AB].*?\(([A-Za-z /[0-9]*?(\(H[0-9]+N[0-9]+\))?)\)')
    FLU_SEG_NUM_REGEX = re.compile(r'Influenza [AB].+ segment ([1-8])')
    KRAKEN_TAX_ID_REGEX = re.compile(r'kraken:taxid\|([0-9]+)\|')
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

        
        

    # @cached_property
    # def data( self ):
    #     """
    #     Returns a dictionary of data from FASTA headers 
        
        
        
    #     Returns:
    #         dictionary with following data:
    #             { 'name':
    #                 { 'segment_number': 
    #                     {
    #                         'fasta_head': 'i|59896322|gb|CY000043|Influenza A......',
    #                         'seq_len': 1234
    #                     }
    #                 }
    #             }
                
    #         Example: if the FASTA header is 
    #         >gi|59896322|gb|CY000043|Influenza A virus (A/New York/34/2003(H3N2)) segment 6, partial sequence
    #         the 'name' of the genome will be: 'Influenza A virus (A/New York/34/2003(H3N2))'
        
    #     """
    #     logging.info( f'scanning file { self.fasta_file_path } for FASTA headers with accepted pattern')
    #     data = {}
    #     with open( self.fasta_file_path ) as fh:
    #         for record in SeqIO.parse(fh, "fasta"):
    #             if not re.search('partial', record.description):
    #                 try:
    #                     match = self.HEADER_REGEX.search( record.description)
    #                     name = match.group(1)
    #                     segment_num = match.group(2)
    #                     seq_len = len(record.seq)
    #                     if name in data:
    #                         if segment_num in data[name]:
    #                             if self.discard_duplicates:
    #                                 continue
    #                             else:
    #                                 raise ValueError(f"found a duplicate definition for name '{name}' segment {segment_num} in {self.fasta_file_path}")
    #                         else:
    #                             data[name][segment_num] = {
    #                                 'fasta_head': record.description,
    #                                 'seq_len': seq_len
    #                             }
    #                     else:
    #                         data[name] = {
    #                             segment_num: {
    #                                 'fasta_head': record.description,
    #                                 'seq_len': seq_len
    #                             }
    #                         }
    #                 except AttributeError:
    #                     # can't parse the header but that is ok, obviously not one we
    #                     # care about, just ignore
    #                     pass
        
    #     logging.info( f'found {len(data.keys())} genomes (unique names) before filtering')
    #     return data