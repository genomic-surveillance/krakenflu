from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import os.path
from cached_property import cached_property
from dataclasses import dataclass
import logging 

from kraken_flu.src.utils import KRAKEN_TAX_ID_REGEX, NCBI_ACC_REGEX
from kraken_flu.src.utils import parse_flu

logging.basicConfig( format='%(asctime)s %(message)s', level=logging.DEBUG )

class FastaHandler():   
    """
    This class handles the genome FASTA files.    
    It can also write the data back to a FASTA file again.
    
    Parameters:
        fasta_file_path: str, required
            path to the FASTA file that is to be parsed
        
    """

    def __init__( self, fasta_file_path: str ):
        if not os.path.isfile( fasta_file_path ):
            raise ValueError(f'file { fasta_file_path } does not exist or is not a file')
        self.fasta_file_path = fasta_file_path

    def _parse_header( self, header:str):
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
            
        flu_type, isolate_name, h_subtype, n_subtype, segment_number = parse_flu( header )
        is_fluA = False
        if flu_type:
            is_flu = True
            if flu_type == 'A':
                is_fluA = True
        else:
            is_flu = False
            
        return flu_type, ncbi_acc, kraken_taxid, is_flu, is_fluA, isolate_name, segment_number, h_subtype, n_subtype

    @cached_property
    def data( self ):
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
            - sequence: DNA sequence
            The following are only applied to flu sequences (None otherwise):
            - is_flu: True if this is a flu sequence
            - flu_name: flu isolate name such as 'A/New York/32/2003(H3N2)' or 'B/Texas/24/2020'
            - flu_seg_num: segment number
            - is_fluA: True if this is an influenza A isolate
            - fluA_H_subtype: H subtype if this is an influenza A isolate, value 'H1', 'H2', etc
            - fluA_N_subtype: N subtype if influenza A, value 'N1', N2', etc.
        
        Returns:
            list of dicts, see above
        
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
                seqlen = len( sequence )
                flu_type, ncbi_acc, kraken_taxid, is_flu, is_fluA, isolate_name, segment_number, h_subtype, n_subtype = self._parse_header(orig_header)
                
                n_all+=1
                if is_flu:
                    if isolate_name is None:
                        unnamed_flu+=1
                    if isolate_name is None:
                        flu_wo_seg_num+=1
                    n_flu+=1

                data.append(
                    FastaRecord(
                        orig_head= orig_header,
                        seqlen= seqlen,
                        sequence = sequence,
                        ncbi_acc= ncbi_acc,
                        is_flu= is_flu,
                        taxid= kraken_taxid,
                        flu_name= isolate_name,
                        flu_seg_num = segment_number,
                        is_fluA= is_fluA,
                        flu_type = flu_type,
                        fluA_H_subtype= h_subtype,
                        fluA_N_subtype= n_subtype) )
        
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
                    fixed_unnamed+=1
                    
            logging.info(f'fixed missing isolate name in {fixed_unnamed} sequences')
        
        
        return data
    
    def write_fasta( self, path:str ):
        """
        Writes FASTA to a new file.
        The output header is identical to the original one, except for influenza sequences, which are written with a 
        modified FASTA header with the format:
        >kraken:taxid|INT|NCBIACC Influenza TYPE (ISOLATE NAME) SEGMENT_NUMBER
        
        Paremeters:
            path: str, required
                Path to the file we are creating
                
        Returns:
            True on success
            
        Side-effects:
            Writes to file
            
        """
        with open( path, 'w' ) as out_fh:
            for record in self.data:
                if record.is_flu:
                    header = ''
                    if record.taxid:
                        header = header + 'kraken:taxid|' + str( record.taxid ) + '|'
                    if record.ncbi_acc:
                        # TODO: not sure if this is needed for kraken2 to assign taxonomy ID later
                        # so just incase, putting the "gb|" back into the NCBI accession ID but only if
                        # it isn't a refseq accession. RefSeq IDs are NOT Genbank IDs, so should not prefix with gb|
                        # https://community.gep.wustl.edu/repository/course_materials_WU/annotation/Genbank_Accessions.pdf
                        if '_'in record.ncbi_acc:
                            ncbi_acc_str = record.ncbi_acc
                        else:
                            ncbi_acc_str = 'gb|' + record.ncbi_acc
                        header = header + ncbi_acc_str
                    header = header + 'Influenza ' + record.flu_type
                    if record.flu_name:
                        header = header + ' (' + record.flu_name + ')'
                    if record.flu_seg_num:
                        header = header + ' segment ' + str( record.flu_seg_num )
                else:
                    header = record.orig_head
                sr = SeqRecord(
                    Seq(record.sequence),
                    id= header,
                    description= ''
                )
                SeqIO.write( sr, out_fh, "fasta" )
        return True
    
@dataclass
class FastaRecord():
    """
    A simple data class to hold the data for a single FASTA record    
    """
    orig_head: str
    seqlen: int
    sequence: str
    ncbi_acc: str
    is_flu: bool
    taxid: int
    flu_name: str
    flu_seg_num: int
    is_fluA: bool
    flu_type: str
    fluA_H_subtype: str
    fluA_N_subtype: str