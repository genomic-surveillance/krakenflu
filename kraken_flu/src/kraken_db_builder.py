import os.path
import logging
from tempfile import NamedTemporaryFile
from collections import defaultdict

from kraken_flu.src.fasta_loader import load_fasta
from kraken_flu.src.taxonomy_loader import load_taxonomy
from kraken_flu.src.db import Db

logging.basicConfig( format='%(asctime)s %(message)s', level=logging.DEBUG )

class KrakenDbBuilder():   
    """
    This class orchestrates the process of creating the modified files for the building of a Kraken DB. 
    It uses the Db class to create a backend database and calls on taxonomy_loader and fasta_loader to 
    populate the DB, then performs the modifications on the taxonomy structure in the DB before exporting the 
    files required to build the final kraken2 DB using the kraken2 build tool.   
    
    TODO: need to check that this example is still correct:
    Consider this example
    of a library build on NCBI viral refseq:
    
    The kraken2 DB preparation commands are:
    kraken2-build \
        --download-taxonomy \
        --db SOME/PATH/
    
    kraken2-build \
        --download-library viral \
        --db SOME/PATH/
        
    This creates the following files in SOME/PATH
    ├── library
    │   └── viral
    │       ├── assembly_summary.txt
    │       ├── library.fna
    │       ├── library.fna.masked
    │       ├── manifest.txt
    │       └── prelim_map.txt
    └── taxonomy
        ├── accmap.dlflag
        ├── citations.dmp
        ├── delnodes.dmp
        ├── division.dmp
        ├── gc.prt
        ├── gencode.dmp
        ├── images.dmp
        ├── merged.dmp
        ├── names.dmp
        ├── nodes.dmp
        ├── nucl_gb.accession2taxid
        ├── nucl_wgs.accession2taxid
        ├── readme.txt
        ├── taxdump.dlflag
        ├── taxdump.tar.gz
        └── taxdump.untarflag
        
    The two required attributes are the full path to 'library' and 'taxonomy'. 
    The files that are read and manipulated are:
    library/library.fna: sequence data
    library/prelim_map.txt: mapping of sequence IDs to taxon IDs
    taxonomy/names.dmp: lists all taxonomy names, will have new names added for flu segment
    taxonomy/nodes.dmp: lists all nodes in the taxonomy, will have new nodes added for flu segments.
    
    When not using pre-built kraken2 reference files, a mapping from NCBI accession ID to taxonomy ID is
    required. This can be obtained from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
    
    Parameters:
        db_path: str, optional
            If provided, the databses will be built at this path and it will not be deleted at the end of 
            the process. If not provided, the DB will be built as a temp file and will be deleted in the end.
        
        db: Db object. optional
            This is only provided for testing and debugging. It allows to pass in a fully built db as a Db object.  
            
    """
    
    def __init__( self, db_path:str=None, db:Db=None):
        self.tax_ids_updated = False
        self.taxonomy_loaded = False
        self.fasta_files_loaded = []
        self.db_path = None
        self._next_new_tax_id = None
        self._new_flu_node_ids = None
        
        if db:
            self._db = db
        else:
            if db_path:
                self.db_path = db_path
            else:
                self.db_path = NamedTemporaryFile().name
            self._db = Db(db_path=self.db_path)
    
    def db_ready(self):
        """
        Returns True if the DB preparations have been completed and the DB is ready for exports
        The DB is ready if we have loaded the taxonomy data and at least one FASTA library. The 
        class has no knowledge of how many FASTA files the user wants to add to the library so it 
        cannot know whether all intended FASTA files have been loaded. 
        """
        if len(self.fasta_files_loaded)>0 and self.taxonomy_loaded:
            return True
        else:
            return False
        
    def next_new_tax_id(self):
        """
        Returns the next available new taxon ID to be used for creating new nodes.  
        Starts off as maximum tax_id in the DB +1, then increments by one every time it is requested
        """
        if not self._next_new_tax_id:
            self._next_new_tax_id = self._db.max_tax_id() + 1
        else:
            self._next_new_tax_id+=1
        
        return self._next_new_tax_id
        
    def load_taxonomy_files(self, taxonomy_dir:str):
        """
        Load taxonomy files from a directory that contains the NCBI taxonomy download.  
        The directory needs to contain files names.dmp and nodes.dmp and may also contain a file 
        of 

        Args:
            taxonomy_dir: str, required
                Path to the directory of downloaded NCBI taxonomy files. Must contain files
                * names.dmp 
                * nodes.dmp  
                If it contains a file "nucl_gb.accession2taxid" then this is also loaded into the DB.
        """
        if not os.path.exists( taxonomy_dir) or not os.path.isdir( taxonomy_dir ):
            raise ValueError(f'path does not exist or is not a directory: {taxonomy_dir}')
        
        names_file_path = os.path.join(taxonomy_dir,'names.dmp')
        if not os.path.exists(names_file_path) or not os.path.isfile(names_file_path):
            raise ValueError(f"cannot find file 'names.dmp' in taxonomy directory {taxonomy_dir}")
        
        nodes_file_path = os.path.join(taxonomy_dir,'nodes.dmp')
        if not os.path.exists(nodes_file_path) or not os.path.isfile(nodes_file_path):
            raise ValueError(f"cannot find file 'nodes.dmp' in taxonomy directory {taxonomy_dir}")
        
        acc2tax_file_path = os.path.join(taxonomy_dir,'nucl_gb.accession2taxid')
        if not os.path.exists( acc2tax_file_path ) or not os.path.isfile( acc2tax_file_path ):
            acc2tax_file_path = None
            
        load_taxonomy(
            db= self._db,  
            names_file_path= names_file_path, 
            nodes_file_path= nodes_file_path, 
            acc2taxid_file_path= acc2tax_file_path
        )
        
        self.taxonomy_loaded= True        

    def load_fasta_file(self, file_path:str, category:str=None):
        """
        Uses the fasta_loader to load a FASTA file into the DB. For details, see fasta_loader module.

        Args:
            file_path: str, required
                Path to the FASTA file we are loading into the DB
            
            category: str, optional
                If provided, will be used to set the column "category" in the "sequences" 
                table, which can be used later to create associations with taxonomy nodes.  
                This is used for cases where we load a specific FASTA file for a known virus 
                (type) and we want to save a hint in the DB for the taxonomy association later. 
        """
        load_fasta(db=self._db, file_path=file_path, category=category)
        self.fasta_files_loaded.append(file_path)
        
    def filter_unnamed_unsegmented_flu(self):
        """
        Marks all flu records as include=0 where the sequence name indicates flu but the isolate  
        name couldn't be parsed, ie the flu_name field is empty and/or segment is empty
        
        Args:
            None
            
        Returns:
            True on success
            
        Side effects:
            sets sequences.include field
        """
        unnamed_flu_sequence_ids = self._db.retrieve_unnamed_unsegmented_flu()
        self._db.mark_as_not_included( unnamed_flu_sequence_ids )
        return True
        
    def filter_incomplete_flu( self, filter_except_patterns:list = [] ):
        """
        Filter out flu genomes that do not have full-length sequences for all 8 segments.
        Genomes that are to be removed are marked by setting attribute include_in_output to False but are 
        not removed from the DB.
        
        Parameters:        
                    
            filter_except_patterns: list(str), optional
                A list of strings that are used to exclude genomes from the Influenza "complete genome" filter.
                Any genome where the FASTA header contains one of the strings in this list will not be subjected to 
                the filter. This was required to ensure that the Goose Guandong H5N1 reference genome (which does not
                have sequences for all 8 segments) is not filtered out.
                
        Returns: 
            number of removed genomes
        
        Side effects: 
            modifies the DB, sets field "include" to False for sequences that are filtered out
        """
        # influenza A/B segment lengths
        # TODO: check that this also applies to influenza C and D, bearing in mind they don't need to 
        # be exact as long as we use the MIN_SEQ_LEN_PROPORTION not too strictly
        MIN_SEG_LENGTHS = {
            1: 2341,
            2: 2341,
            3: 2233,
            4: 1778,
            5: 1565,
            6: 1413,
            7: 1027,
            8: 890
        }
        
        # minimum proportion of each of the segments that must be covered
        MIN_SEQ_LEN_PROPORTION = 0.9
        
        logging.info( f'starting filter to remove incomplete flu genomes')
    
        flu_data = self._db.get_flu_name_segment_data_dict()
        sequence_ids_to_remove = []
        isolate_removal_count = 0
        isolate_total_count = 0
        isolate_exempt_count = 0
        for isolate_name, segment_data in flu_data.items():
            isolate_total_count+=1
            if isolate_name in filter_except_patterns:
                # this one should not be subjected to filtering - skip
                isolate_exempt_count+=1
                continue
            remove_isolate=False
            # filter out if we don't have all 8 segments
            if sorted(segment_data.keys()) != list(range(1,9)):
                remove_isolate = True
            for segment_number, length in segment_data.items():
                if length < MIN_SEG_LENGTHS[ segment_number ] * MIN_SEQ_LEN_PROPORTION:
                    remove_isolate = True
                    continue
            if remove_isolate:
                ids = self._db.retrieve_sequence_ids_by_flu_name(isolate_name)
                sequence_ids_to_remove.extend(ids)
                isolate_removal_count+=1
        
        logging.info( f'{isolate_removal_count} of {isolate_total_count} flu isolates were identified as incomplete ({isolate_exempt_count} matched exempt list). Marking in database')
        self._db.mark_as_not_included(sequence_ids_to_remove)
        logging.info( f'Completed updating database')
        return True
    
    def create_segmented_flu_taxonomy_nodes(self, missing_parent_is_exception:bool=False):
        """
        Creates the new taxon nodes for flu and returns tax_ids datastructure.  
        The result is a two new levels of "artificial" taxon nodes that looks like this for influenza A: 
                                    ┌──────────────────┐                                                         
                                    │ Influenza A virus│                                                         
                   ┌────────────────┴──────┬───────────┴──────────────────┐                                  
        ┌──────────▼──────────┐ ┌──────────▼──────────┐       ┌───────────▼──────────┐                       
        │Influenza A segment 1│ │Influenza A segment x│ [...] │Influenza A segment 4 │ [...]                 
        └─┬───────────────────┘ └─────────────────────┘  ┌────┴──────────────────┬───┘                        
          │ ┌──────────────────────────┐   ┌─────────────▼──────────┐ ┌──────────▼─────────────┐             
          ├─► isolate segment sequence │   │Influenza A H1 segment 4│ │Influenza A H2 segment 4│ [...]       
          │ └──────────────────────────┘   └┬───────────────────────┘ └──────┬─────────────────┘             
          │ ┌──────────────────────────┐    │  ┌──────────────────────────┐  │  ┌──────────────────────────┐ 
          └─► isolate segment sequence │    ├──► isolate segment sequence │  ├──► isolate segment sequence │ 
            └──────────────────────────┘    │  └──────────────────────────┘  │  └──────────────────────────┘ 
                [...]                       │  ┌──────────────────────────┐  │  ┌──────────────────────────┐ 
                                            └──► isolate segment sequence │  └──► isolate segment sequence │ 
                                               └──────────────────────────┘     └──────────────────────────┘ 
                                                    [...]                               [...]              
        For flu B,C and D, the additional level for Hx and Nx segments are not created because these flu types 
        are not further classified by segments 4 and 6.  
        
        The creation of new nodes will only be run once. If it has already been run, the previously generated 
        new_node_ids data is returned as-is.  

        Parameters:
            missing_parent_is_exception: bool, optional, defaults to False
                Parental nodes for Influenza X Virus must exist in the DB for all virus types prior to 
                building the new taxonomy structure.  
                If this flag is True, a missing type is treated as an exception, otherwise it is simply 
                shown in the logs and the building ot the custom taxonomy is skipped for that virus type
        
        Returns:
            dict of new node names to node taxon ids
            Structure:
            { tpye: { subtype: {segment_number: tax_id}}}
        
        """
        if self._new_flu_node_ids:
            return self._new_flu_node_ids
        
        types = ['A','B','C','D']
        new_node_ids = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        for type in types:
            parent_node_name = f"Influenza {type} virus"
            parent_tax_id = self._db.retrieve_tax_id_by_node_scientific_name(parent_node_name)
            if not parent_tax_id:
                if missing_parent_is_exception:
                    raise ValueError(f"parent node '{parent_node_name}' not found in DB")
                else:
                    logging.warning(f"parent node '{parent_node_name}' not found in DB - skipping")
                    continue
            # create the new type/segment nodes
            for seg_num in range(1,9):
                new_tax_id = self.next_new_tax_id()
                subtype_str= None
                self._db.add_taxon( tax_id= new_tax_id, parent_tax_id= parent_tax_id, name=f'Influenza {type} segment {seg_num}' )
                new_node_ids[type][subtype_str][seg_num] = new_tax_id

                # create all Hx and Nx segment nodes if this is flu A and one of the relevant segments (4 or 6)
                if type=='A':
                    # the nodes created here should have the segment node as its parent
                    segment_parent_tax_id = new_tax_id
                    if seg_num == 4:
                        subtype_letter = 'H'
                        max_subtype_num = 18
                    elif seg_num == 6:
                        subtype_letter = 'N'
                        max_subtype_num = 11
                    else:
                        continue
                    for subtype_num in range(1,max_subtype_num+1):
                        new_tax_id = self.next_new_tax_id()
                        subtype_str = subtype_letter+str(subtype_num)
                        self._db.add_taxon( tax_id= new_tax_id, parent_tax_id= segment_parent_tax_id, name=f'Influenza {type} {subtype_str} segment {seg_num}' )
                        new_node_ids[type][subtype_str][seg_num] = new_tax_id
        self._new_flu_node_ids = new_node_ids
        return new_node_ids        
                
        
    def update_fasta_tax_ids( self ):
        """
        Combines the TaxonomyHandler and FastaHandler to update taxon IDs for output of the new
        FASTA file.
        
        Paramters:
            none
            
        Returns:
            True on success
            
        Side effects:
            - updates taxids in slef._fasta_handler.data
            - sets self.tax_ids_update to True
        """
        if self.tax_ids_updated:
            # already done this, nothing to do
            return True
                
        logging.info( f'starting to assign new taxonomy IDs to sequence records')
        
        # trigger the cascade to create the new taxa
        self._taxonomy_handler.create_influenza_isolate_segment_taxa()
        
        for fasta_record in self._fasta_handler.data:
            # currently, we only make changes to flu sequences
            if fasta_record.is_flu and fasta_record.flu_name and fasta_record.flu_seg_num:
                try:
                    taxid = self._taxonomy_handler.influenza_isolate_segment_tax_ids[ fasta_record.flu_name ][ fasta_record.flu_seg_num ]
                    fasta_record.taxid = taxid
                except KeyError:
                    logging.info( f'found influenza record for isolate "{fasta_record.flu_name}" segment { fasta_record.flu_seg_num} in FASTA file but no corresponding entry in taxonomy - skipping')

        logging.info( f'finished assigning new taxonomy IDs')

        self.tax_ids_updated = True
        return True

    def create_db_ready_dir( self, path: str, force: bool = True, fasta_file_name:str = 'library.fna', filter_incomplete_flu:bool = True,  filter_except_patterns: list= [], drop_unparsed_flu:bool = True ):
        """
        Create a directory with all files needed to build a new kraken database. 
        This method triggers the data acquisition from the FastaHandler and TaxonomyHandler, then runs the
        changing of taxonomy IDs (delegating to the two handlers) before writing a new FASTA file with updated 
        kraken:taxid tags as well as the names and nodes files of the NCBI taxonomy with the added taxa.
        
        The FASTA file is named library.fna as per kraken2 convention but can be given a different name. The 
        NCBI taxonomy files names.dmp and nodes.dmp cannot be renamed as kraken2 relies on those names.
        
        Parameters:
            path: str
                path to the directory into which the files will be written.
                If force is not used, the DIR must not exist yet
                
            force: bool, optional, defaults to True
                if true, force overwrite whatever is already in the output dir. If not, will fail
                if directory already exists
                
            filter_incomplete_flu: bool, optional, defaults to True
                If True, runs the filter that removes incomplete influenza genomes. See FastaHandler for details
                            
            filter_except_patterns: list(str), optional
                A list of strings that are used to exclude genomes from the Influenza "complete genome" filter.
                Any genome where the FASTA header contains one of the strings in this list will not be subjected to 
                the filter. This was required to ensure that the Goose Guandong H5N1 reference genome (which does not
                have sequences for all 8 segments) is not filtered out.
                
            drop_unparsed_flu: bool, optional, defaults to True
                If True, sequences whose name indicates that it is flu but where the name cannot be parsed
                properly (we can't obtain an isolate name) will be dropped. Without this filter, such cases
                would end up as flu whole genomes that are not segmented like the rest of the flu genomes.
            
        Returns:
            True if success
            
        Side effects:
            Writes files to path
            
        """
        library_path = os.path.join( path , 'library' )
        taxonomy_path = os.path.join( path, 'taxonomy' )
        
        fasta_file_path = os.path.join( library_path, fasta_file_name )
        names_file_path = os.path.join( taxonomy_path, 'names.dmp')
        nodes_file_path = os.path.join( taxonomy_path, 'nodes.dmp')

        if os.path.exists( path ): 
            if not force:
                raise ValueError(f'directory { path } exists already. Will not write into existing directory')
        else:
            os.mkdir( path )
            
        if not os.path.exists( library_path ):
            os.mkdir( library_path )
        if not os.path.exists( taxonomy_path ):
            os.mkdir( taxonomy_path )
        
        logging.info( f'found { self._fasta_handler.n_seq_total()} sequence records in {self.fasta_file_path}')
        logging.info( f'{ self._fasta_handler.n_seq_flu()} sequence records identified as influenza')
        
        if filter_incomplete_flu or drop_unparsed_flu:
            logging.info( f'starting to apply filters on influenza genomes')
            n_seq_filtered = self._fasta_handler.n_seq_filtered() # should be 0 at this point
            
            if filter_incomplete_flu:
                logging.info( f'starting to filter incomplete influenza genomes')
                self._fasta_handler.remove_incomplete_flu( filter_except_patterns= filter_except_patterns )
                n_seq_filtered = self._fasta_handler.n_seq_filtered() - n_seq_filtered
                logging.info( f'{ n_seq_filtered} flu sequence records have been marked to be removed from output due to being incomplete')
            if drop_unparsed_flu:
                logging.info( f'starting to apply drop_unparsed_flu filter')
                self._fasta_handler.remove_unparsed_flu()
                n_seq_filtered = self._fasta_handler.n_seq_filtered() - n_seq_filtered
                logging.info( f'{ n_seq_filtered } flu sequence records have been marked to be removed from output as unparsable')
        
            logging.info( f'after applying filters, { self._fasta_handler.n_seq_filtered()} flu sequence records have been marked to be removed from output')
        
        if not self.tax_ids_updated:
            self.update_fasta_tax_ids()
        
        logging.info( f'writing output files')
        self._fasta_handler.write_fasta( fasta_file_path )
        self._taxonomy_handler.write_names_file( names_file_path )
        self._taxonomy_handler.write_nodes_file( nodes_file_path )
        logging.info( f'---- process complete ----')
        return True

