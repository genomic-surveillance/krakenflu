from Bio import SeqIO
import re
import os.path
from cached_property import cached_property
import logging 

logging.basicConfig( format='%(asctime)s %(message)s', level=logging.DEBUG )

class NcbiInfluenzaFastaFilter():   
    """
    This class handles the filtering of large-scale downloads from 
    https://ftp.ncbi.nih.gov/genomes/INFLUENZA/influenza.fna
    which contains all influenza data from NCBI.
    
    Filters the file down the sequences tht fulfil the following criteria:
    
        - genome has all 8 segments, each with at least 90% of the expected sequence length
        - segments must be named according to convention such as this example
            "Influenza A virus (A/ruddy turnstone/Delaware Bay/205/2017(H3N2)) segment 1"
    
    Parameters:
        fasta_file_path: str, required
            path to the FASTA file that is to be filtered
            
        discard_duplicates: bool, optional, defaults to False
            if True, discard FASTA headers where name and segment number is duplicated in the file.
            This happens for partial sequences where separate sequence records are provided for the
            different genes on the same segment. Such sequences are not of interest for us as they are
            not complete genome segments anyway
        
    """
    
    HEADER_REGEX = re.compile(r'(Influenza [AB].+) segment ([1-8])')
    
    # influenza A/B segment lengths
    MIN_SEG_LENGTHS = {
        '1': 2341,
        '2': 2341,
        '3': 2233,
        '4': 1778,
        '5': 1565,
        '6': 1413,
        '7': 1027,
        '8': 890
    }
    
    # minimum proportion of each of the segments that must be covered
    MIN_SEQ_LEN_PROPORTION = 0.9
    
    def __init__( self, fasta_file_path: str, discard_duplicates: bool = False ):
        if not os.path.isfile( fasta_file_path ):
            raise ValueError(f'file { fasta_file_path } does not exist or is not a file')
        self.fasta_file_path = fasta_file_path
        self.discard_duplicates = discard_duplicates
        
    @cached_property
    def records_by_name( self ):
        """
        Returns a dictionary of data from FASTA headers with the name of the virus used as the key, 
        thus grouping the sequences into genomes.
        NOTE: this relies on unique names being given by source labs. This SHOULD be the case if the
        source lab adheres to the conventional influenza naming. 
        If needed, a more robust approach would be to require the NCBI acc2taxid file, look up each gb ID
        in that file and use the tax ID to group sequences into genomes.
        
        Returns:
            dictionary with following data:
                { 'name':
                    { 'segment_number': 
                        {
                            'fasta_head': 'i|59896322|gb|CY000043|Influenza A......',
                            'seq_len': 1234
                        }
                    }
                }
                
            Example: if the FASTA header is 
            >gi|59896322|gb|CY000043|Influenza A virus (A/New York/34/2003(H3N2)) segment 6, partial sequence
            the 'name' of the genome will be: 'Influenza A virus (A/New York/34/2003(H3N2))'
        
        """
        logging.info( f'scanning file { self.fasta_file_path } for FASTA headers with accepted pattern')
        data = {}
        with open( self.fasta_file_path ) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                if not re.search('partial', record.description):
                    try:
                        match = self.HEADER_REGEX.search( record.description)
                        name = match.group(1)
                        segment_num = match.group(2)
                        seq_len = len(record.seq)
                        if name in data:
                            if segment_num in data[name]:
                                if self.discard_duplicates:
                                    continue
                                else:
                                    raise ValueError(f"found a duplicate definition for name '{name}' segment {segment_num} in {self.fasta_file_path}")
                            else:
                                data[name][segment_num] = {
                                    'fasta_head': record.description,
                                    'seq_len': seq_len
                                }
                        else:
                            data[name] = {
                                segment_num: {
                                    'fasta_head': record.description,
                                    'seq_len': seq_len
                                }
                            }
                    except AttributeError:
                        # can't parse the header but that is ok, obviously not one we
                        # care about, just ignore
                        pass
        
        logging.info( f'found {len(data.keys())} genomes (unique names) before filtering')
        return data
        
    @cached_property
    def filtered_fasta_headers(self):
        """
        Applies the filters to the extracted records and returns a list of the FASTA headers that
        should be kept, where the header is the record.description from SeqIO
        
        Returns:
            list of FASTA headers to keep
        """
        to_keep = []
        for name, data in self.records_by_name.items():
            if [data.keys()].sort() == [1, 2, 3, 4, 5, 6, 7, 8].sort():
                # we have all 8 segments, check if any of them are too short
                full_length_headers = []
                for segment in data.keys():
                    if data[segment]['seq_len'] >= self.MIN_SEG_LENGTHS[segment] * self.MIN_SEQ_LEN_PROPORTION:
                        full_length_headers.append( data[segment]['fasta_head'])
                
                if len( full_length_headers ) == 8:
                    to_keep.extend (full_length_headers)  
        
        logging.info(f'{len(to_keep)} FASTA records pass the filters and will be written to output file')
        
        return to_keep
        
    def write_filtered_file( self, out_path: str ):
        """
        Writes the final, filtered, file to out_path
        
        Parameters:
            out_path: str, required
                path to the file which will be created
        """
        dir = str(os.path.dirname( out_path ))
        if not os.path.isdir( dir ):
            raise ValueError(f'directory "{dir}" does not exist, cannot create file "{out_path}"')
        
        with open( out_path, 'w') as out_fh:
            with open( self.fasta_file_path ) as in_fh:
                for record in SeqIO.parse(in_fh, "fasta"):
                    if record.description in self.filtered_fasta_headers:
                        out_fh.write( '>' + record.description +"\n")
                        out_fh.write( str(record.seq) + "\n")
        logging.info( f'finished writing filtered genomes to {out_path}')