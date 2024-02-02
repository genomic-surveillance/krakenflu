from Bio import SeqIO
import re
import os.path
from cached_property import cached_property
from dataclasses import dataclass
import logging 

from kraken_flu.src.utils import FLU_REGEX,FLU_ISOLATE_NAME_REGEX, FLU_SEG_NUM_REGEX, KRAKEN_TAX_ID_REGEX, NCBI_ACC_REGEX

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
            path to the FASTA file that is to be parsed
        
    """

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
            ncbi_acc = NCBI_ACC_REGEX.search( header ).group(0)
            ncbi_acc = re.sub(r'^gb\|', '', ncbi_acc)
        except AttributeError:
            # an NCBI accession ID is required
            raise ValueError(f"could not parse NCBI acc ID from FASTA header '{header}'")

        if (match := KRAKEN_TAX_ID_REGEX.search( header )) is not None:
            kraken_taxid = int(match.group(1))
        else:
            kraken_taxid = None
            
        flu_isolate_name = None
        flu_segment_number = None
        if FLU_REGEX.search( header ):
            is_flu = True
            if (match := FLU_SEG_NUM_REGEX.search( header )) is not None:
                flu_segment_number = int(match.group(1))
            if (match := FLU_ISOLATE_NAME_REGEX.search( header )) is not None:
                flu_isolate_name = match.group(1)
        else:
            is_flu = False
            
        return ncbi_acc, kraken_taxid, is_flu, flu_isolate_name, flu_segment_number 

    @cached_property
    def data( self ):
        """
        Returns FASTA data as a list of dicts, each dict has the following keys:
            - orig_head: the original un-modified FASTA header
            - sequence: the DNA sequence
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
            
        TODO:
        Reading the sequence data into RAM because it is easier to work with it this way and even 
        the largest file we are currently working with (all influenza from NCBI) is only around 500MB, 
        which can easily be handled on the farm. But this might become a problem if we need to work with 
        larger files and might need to change.
        
        """
        logging.info( f'parsing FASTA headers from { self.fasta_file_path }')
        data = []
        n_all = 0
        n_flu = 0
        unnamed_flu = 0
        flu_wo_seg_num = 0
        with open( self.fasta_file_path ) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                orig_header = record.description
                sequence = record.seq
                seqlen = len(sequence)
                ncbi_acc, kraken_taxid, is_flu, flu_isolate_name, flu_segment_number = self._parse_header(orig_header)
                
                # TODO: not sure if this is needed for kraken2 to assign taxonomy ID later
                # so just incase, putting the "gb|" back into the NCBI accession ID but only if
                # it isn't a refseq accession. RefSeq IDs are NOT Genbank IDs, so should not prefix with gb|
                # https://community.gep.wustl.edu/repository/course_materials_WU/annotation/Genbank_Accessions.pdf
                if '_'in ncbi_acc:
                    ncbi_acc_str = ncbi_acc
                else:
                    ncbi_acc_str = 'gb|'+ncbi_acc
                    
                n_all+=1
                if is_flu:
                    if flu_isolate_name is None:
                        unnamed_flu+=1
                    if flu_segment_number is None:
                        flu_wo_seg_num+=1
                    n_flu+=1
                    mod_header = ' '.join([
                        ncbi_acc_str +'|', 
                        'Influenza', 
                        flu_isolate_name or 'unnamed', 
                        'segment', 
                        str(flu_segment_number)])
                else:
                    mod_header = orig_header

                data.append(
                    FastaRecord(
                        orig_head= orig_header,
                        sequence= sequence,
                        mod_head= mod_header,
                        seqlen= seqlen,
                        ncbi_acc= ncbi_acc,
                        is_flu= is_flu,
                        taxid= kraken_taxid,
                        flu_name= flu_isolate_name,
                        flu_seg_num = flu_segment_number ) )
        
        logging.info( f'found {n_all} sequences, {n_flu} of which are influenza')
        if flu_wo_seg_num > 0:
            logging.info(f'found {flu_wo_seg_num} influenza sequences without a segment number')
            
        if unnamed_flu>0:
            logging.info(f'found {unnamed_flu} influenza sequences without an isolate name - attempting to derive from taxid')
            
            # this fix attempts to fill in missing isolate names fromm other sequences that
            # have the same tax id. This fixes the problem with the Influenza B RefSeq, which 
            # has segment 1 annotated with a non-standard FASTA header that is missing the isolate name
            taxid2isolate_name= {}
            for d in data:
                if d.flu_name is not None and d.taxid is not None:
                    taxid2isolate_name[ d.taxid ] = d.flu_name 
                
            fixed_unnamed= 0
            for d in data:
                if d.is_flu and d.taxid is not None and d.flu_name is None and d.taxid in taxid2isolate_name:
                    d.flu_name = taxid2isolate_name[ d.taxid ]
                    d.mod_head = re.sub('unnamed', d.flu_name, d.mod_head)
                    fixed_unnamed+=1
                    
            logging.info(f'fixed missing isolate name in {fixed_unnamed} sequences')
        
        
        return data
    
@dataclass
class FastaRecord():
    """
    A simple data class to hold the data for a single FASTA record    
    """
    orig_head: str
    sequence: str
    mod_head: str
    seqlen: int
    ncbi_acc: str
    is_flu: bool
    taxid: int
    flu_name: str
    flu_seg_num: int
        
    
    