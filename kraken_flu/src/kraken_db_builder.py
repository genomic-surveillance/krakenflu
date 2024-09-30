import os.path
import logging
from tempfile import NamedTemporaryFile
from collections import defaultdict
from math import floor

from kraken_flu.src.fasta_loader import load_fasta
from kraken_flu.src.fasta_writer import write_fasta
from kraken_flu.src.taxonomy_loader import load_taxonomy
from kraken_flu.src.taxonomy_writer import write_taxonomy
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
        
    def load_taxonomy_files(self, taxonomy_dir:str, no_acc2taxid:bool=False):
        """
        Load taxonomy files from a directory that contains the NCBI taxonomy download.  
        The directory needs to contain files names.dmp and nodes.dmp and may also contain a file 
        of accession to taxonomy ID mappings "nucl_gb.accession2taxid"

        Args:
            taxonomy_dir: str, required
                Path to the directory of downloaded NCBI taxonomy files. Must contain files
                * names.dmp 
                * nodes.dmp  
                If it contains a file "nucl_gb.accession2taxid" then this is also loaded into the DB unless option 
                no_acc2taxid is in use.
                
            no_acc2taxid: bool, optional, defaults to False
                if True, ignore any NCBI acc2taxid file in the taxonomy directory.  The file 
                will not be loaded into the kraken_flu DB and no general linking of sequences to 
                taxa will take place, ie only the special cases such as flu will end up hacing kraken:taxid tags 
                in the final FASTA file and all other sequences need to be linked by the kraken2 
                build process.  
                
        """
        if not os.path.exists( taxonomy_dir) or not os.path.isdir( taxonomy_dir ):
            raise ValueError(f'path does not exist or is not a directory: {taxonomy_dir}')
        
        names_file_path = os.path.join(taxonomy_dir,'names.dmp')
        if not os.path.exists(names_file_path) or not os.path.isfile(names_file_path):
            raise ValueError(f"cannot find file 'names.dmp' in taxonomy directory {taxonomy_dir}")
        
        nodes_file_path = os.path.join(taxonomy_dir,'nodes.dmp')
        if not os.path.exists(nodes_file_path) or not os.path.isfile(nodes_file_path):
            raise ValueError(f"cannot find file 'nodes.dmp' in taxonomy directory {taxonomy_dir}")
        
        # The accession-to-taxid file is optional. If it is present, we will load it and use it to link
        # all FASTA records to a taxid (like kraken2 build does). If it is not present or an option was used 
        # to ignore it, we will skip this step. This will mean that only special cases, like flu, will have kraken:taxid 
        # tags in the FASTA header at the end of the process. All other sequences must be linked by kraken2 build.
        acc2tax_file_path = os.path.join(taxonomy_dir,'nucl_gb.accession2taxid')
        if os.path.isfile( acc2tax_file_path ):
            if no_acc2taxid:
                logging.info(f"option no_acc2taxid in use: file {acc2tax_file_path} found in taxonomy dir but will be ignored" )
                acc2tax_file_path = None
            else:
                logging.info(f"found acc2txid file {acc2tax_file_path} in taxonomy dir - will load and use NCBI accession-to-taxid data. Use option --no-acc2taxid to avoid this.")
        else:
            logging.info("could not find an NCBI accession-to-taxid file in taxonomy dir. Skipping the general linking of sequences to taxids. Check manual for further details.")
            acc2tax_file_path= None
                
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
        This method is used for "generic" sequences that do not require a category label. This includes 
        the general RefSeq FASTA as well as influenza data, where the FASTA headers are used to figure out 
        which taxonomy node to link to.  
        For sequences that require a category label (such as RSV A/B), use the load_labelled_fasta method instead.  

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
        logging.info("started excluding unnamed and/or unsegmented influenza records from DB")
        unnamed_flu_sequence_ids = self._db.retrieve_unnamed_unsegmented_flu()
        self._db.mark_as_not_included( unnamed_flu_sequence_ids )
        logging.info("finished excluding unnamed and/or unsegmented influenza records")
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
        n_total = len(flu_data)
        logging.info( f"found {n_total} distinct influenza isolate names in the DB - starting to analyse them for completeness")
        
        sequence_ids_to_remove = []
        isolate_removal_count = 0
        isolate_total_count = 0
        isolate_exempt_count = 0
        last_log_cp = 0
        for isolate_name, segment_data in flu_data.items():
            isolate_total_count+=1
            percent_done = ( isolate_total_count / n_total ) * 100
            log_cp = int(floor(percent_done /10 ))
            if log_cp > last_log_cp:
                last_log_cp = log_cp
                logging.info(f"{log_cp * 10 }% complete")
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
    
    def filter_flu_a_wo_subtype(self):
        """
        Exclude all flu A sequences where the H and N subtypes are not set.  
        If we don't do this, flu A sequences will become linked to "Influenza A segment 4" 
        and "Influenza A segment 6" nodes, thus bypassing the "Influenza A Hx segment 4" and 
        "Influenza A Nx segment 4" level of taxonomy nodes created by create_segmented_flu_taxonomy_nodes.  
        The main source of flu A genomes in this category are those that have a "mixed" subtype in the 
        isolate name.   
        TODO: investigate whether we should create nodes "Infleunza A mixed subtype segment 4" (and 6) 
        and keep the mixed subtype sequences. They might be be valuable for the analysis but this will 
        need some research.  
        
        Args:
            none
            
        Side effects:
            Sets sequences.include field for some flu A sequences
            
        Returns:
            True on success
        """
        logging.info('starting to filter out flu A genomes without subtype (such as "mixed" subtypes)')
        sequence_ids_to_remove = self._db.retrieve_ids_flu_a_wo_subtype()
        self._db.mark_as_not_included(sequence_ids_to_remove)
        logging.info(f'filtered out {len(sequence_ids_to_remove)} flu A genomes without a subtype')
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
        
        logging.info("starting to create new taxonomy node for segmented flu genomes")

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
        logging.info("created new taxonomy node for segmented flu genomes")
        return new_node_ids        
            
    def assign_flu_taxonomy_nodes(self):
        """
        Assigns sequences that are identified as flu to the new taxonomy nodes from method 
        "create_segmented_flu_taxonomy_nodes", which is run at this point if it hasn't been run already.  
        This involves the following steps:
            - create a new taxonomy node for the flue genome in the taxonomy
            - set the parent of this new node to one of the nodes created by  "create_segmented_flu_taxonomy_nodes"
            - set the tax_id of the sequence record to the tax_id of the new record in taxonomy_nodes

        NOTE: The rationale for creating a new taxonomy node for each genome is that we provide multiple 
        genomes for influenza. In order to make kraken2 assign reads to a specific genome, each genome must have its own 
        node in the taxonomy. This is not the case for RefSeq sequences in the NCBI taxonomy, because there is usually 
        only one genome for each species.  
        For these reasons, every flu genome in our modified taxonomy has its own tax_id, whereas the flu A genome segments 
        downloaded from NCBI all have the same tax_id, which is the tax_id of "Influenza A virus"
        
        Returns:
            True on success
            
        Side effect:
            - creates new records in taxonomy_nodes and taxonomy_names
            - sequences.tax_id
        """
        # get the tax_ids for the new nodes or run the creation of the new nodes if not done already
        if not self._new_flu_node_ids:
            new_flu_node_ids = self.create_segmented_flu_taxonomy_nodes()
        else:
            new_flu_node_ids = self._new_flu_node_ids
            
        logging.info("starting to assign flu genomes to new taxonomy nodes for segmented flu")
        
        flu_sequences = self._db.retrieve_all_flu_sequences(skip_excluded= True)
        n_records = len(flu_sequences)
        n=0
        last_log_cp = 0
        logging.info(f"retrieved {n_records} flu sequence records from database")
        
        new_tax_id = self.next_new_tax_id()
        with self._db.bulk_update_buffer(table_name='sequences', id_field='id', update_fields= ['tax_id','mod_fasta_header'], buffer_size= 50000) as seq_update_buffer:

            with self._db.bulk_insert_buffer(table_name='taxonomy_nodes', buffer_size= 50000) as nodes_buffer:
                with self._db.bulk_insert_buffer(table_name='taxonomy_names', buffer_size= 50000) as names_buffer:
                    for flu_sequence in flu_sequences:
                        new_tax_id+=1
                        flu_name= flu_sequence['flu_name']
                        flu_type= flu_sequence['flu_type']
                        segment_number= flu_sequence['segment_number']
                        if not flu_name or not segment_number:
                            continue
                        flu_a_h_subtype= flu_sequence['flu_a_h_subtype']
                        flu_a_n_subtype= flu_sequence['flu_a_n_subtype']
                        id= flu_sequence['id']
                        n+=1
                        percent_done = ( n / n_records ) * 100
                        log_cp = int(floor(percent_done /10 ))
                        if log_cp > last_log_cp:
                            last_log_cp = log_cp
                            logging.info(f"{log_cp * 10 }% complete")
                        if segment_number==4 and flu_a_h_subtype:
                            subtype='H' + str(flu_a_h_subtype)
                        elif segment_number==6 and flu_a_n_subtype:
                            subtype='N' + str(flu_a_n_subtype) 
                        else:
                            subtype= None
                            
                        # provide a unified alternative name for the sequence that will be written to the 
                        # final FASTA output file
                        mod_fasta_header = flu_name + ' segment ' + str(segment_number)
                            
                        # Get the tax_id of the matching parent node created by "create_segmented_flu_taxonomy_nodes"
                        parent_tax_id = new_flu_node_ids[flu_type][subtype][segment_number]
                        if not parent_tax_id:
                            raise ValueError(f"could not find the parent node for flu type: {flu_sequence['flu_type']}, subtype: {subtype}, segment: {flu_sequence['segment_number']}")    
                            
                        name= ' '.join([flu_name, 'segment',str(segment_number)])
                        # create the new node/name records for the flu sequence and link it to the parent node
                        n_inserted_nodes = nodes_buffer.add_row(
                            {
                                'tax_id': new_tax_id,
                                'parent_tax_id':parent_tax_id,
                                'rank': 'no rank',
                                'embl_code': None,
                                'division_id': 9,                   
                                'inherited_div_flag': 1,            
                                'genetic_code_id': 1,
                                'inherited_GC_flag': 1,
                                'mitochondrial_genetic_code_id': 0,
                                'inherited_MGC_flag': 1,
                                'GenBank_hidden_flag': 0,
                                'hidden_subtree_root_flag': 0,
                                'comments': 'kraken_flu added node'
                            }
                        )
                        if n_inserted_nodes > 0:
                            logging.info(f'flushed {n_inserted_nodes} nodes records to DB')

                        n_inserted_names = names_buffer.add_row(
                            {
                                'tax_id': new_tax_id,
                                'name': name,
                                'unique_name': name,
                                'name_class': 'scientific name'
                            }
                        )
                        if n_inserted_names > 0:
                            logging.info(f'flushed {n_inserted_names} names records to DB')
                            
                        # link sequence records to the newly inserted taxonomy nodes
                        n_updated_seqs= seq_update_buffer.add_row(
                            {
                                'id': id,
                                'tax_id': new_tax_id,
                                'mod_fasta_header': mod_fasta_header,
                            }
                        )
                        if n_updated_seqs > 0:
                            logging.info(f'flushed {n_updated_seqs} sequence record updates to DB')
        
        logging.info("finished setting taxonomy IDs for segmented flu genomes")
        return True

    def create_rsv_taxonomy_from_files(self, rsv_size_filter:bool=False):
        """
        Discard RSV sequences except those that were explicitly loaded into the database as RSV A 
        and RSV B genome sequences. These are labelled 'RSV A' and 'RSV B' in the sequences.category field.  
        Only these explicitly loaded RSV A/B type sequences are linked to the respective hRSV A and B parent nodes. 
        The parent nodes are identified by name.  
        
        An optional size filter can be used to remove incomplete genomes.  
        
        Rationale for this method:
        RSV sequences in RefSeq are currently not sufficient for our purposes. For this reason, it was decided 
        that we would obtain RSV sequences classified as A or B from Nextstrain instead. We obtain these as 
        separate downloads, which can be loaded into the kraken-flu database with labels "RSV A" and "RSV B" 
        which are stored in the sequences.category field. The aim of this method is to replace the existing hRSV 
        taxonomy bny removing all sequences that are linked to the hRSV parent taxonomy node, then replacing them 
        with the Nextstrain (or other authority) typed sequences based on the category labels in the database.  
        If we would skip the removal of existing hRSV sequences, we would end up with sequences that are directly 
        linked to the hRSV parent node, which bypasses our A/B type classification.  

        Args:
            rsv_size_filter: bool, defaults to False
                If True, RSV sequences are filtered on size to keep only full-length or nearly 
                full length genomes.  
        """
        for t in (['A','B']):
            label = 'RSV ' + t
            if not self._db.sequences_category_exists(label):
                raise ValueError(f'method create_rsv_taxonomy_from_file called but no sequences loaded into DB with label "{label}"')
        
        if rsv_size_filter:
            self._apply_rsv_size_filter(categories=['RSV A','RSV B'])
        
        hrsv_parent_tax_id = self._db.retrieve_tax_id_by_node_scientific_name('human respiratory syncytial virus')
        if not hrsv_parent_tax_id:
            raise ValueError('could not find taxonomy node for "human respiratory syncytial virus" in database')
        
        hrsv_a_parent_tax_id = self._db.retrieve_tax_id_by_node_scientific_name('Human respiratory syncytial virus A')
        if not hrsv_a_parent_tax_id:
            raise ValueError('could not find taxonomy node for "Human respiratory syncytial virus A" in database')
        
        hrsv_b_parent_tax_id = self._db.retrieve_tax_id_by_node_scientific_name('Human respiratory syncytial virus B')
        if not hrsv_b_parent_tax_id:
            raise ValueError('could not find taxonomy node for "Human respiratory syncytial virus B" in database')
            
        # Remove all sequences that are currently linked to the parent hRSV node and all nodes below it.  
        # We want to replace the existing RSV sequences with the ones labelled specifically as A or B type and
        # if we would leave the existing sequences in the taxonomy tree, kraken2 would map reads to them which would 
        # greatly reduce the number of useful matches that are directly assigned to RSV A or B
        self.filter_out_sequences_linked_to_taxonomy_sub_tree(tax_id= hrsv_parent_tax_id)
            
    def filter_out_sequences_linked_to_taxonomy_sub_tree(self, tax_id:int):
        """
        Filter out (set sequences.include=0) all sequences linked to a taxonomy sub-tree, starting at a
        parent taxonomy node, identified by a tax_id.  
        The method traverses the sub-tree recursively and identifies all sequences that have a tax_id linkage 
        to any of the nodes in the sub-tree. These are filtered by setting the include field.  

        Args:
            tax_id: int, required
                The tax_id of the parent node of the sub-tree from which we want to filter out linked 
                sequences. All sequences linked to this node or any of its children are filtered.  
                
        Returns:
            list of sequences.id for the sequences that were removed
        """
        sequence_ids_to_remove = self._db.get_sequence_ids_linked_to_taxon(tax_id= tax_id, include_children= True)
        self._db.mark_as_not_included(sequence_ids_to_remove)
        return sequence_ids_to_remove
            
    def _apply_rsv_size_filter(self, categories:list=None):
        """
        Filter out (mark sequences.includee=0) all RSV genomes that do not meet the minimum length filter.  
        This currently relies on the sequences.cateogry label being set to 'RSV A' or 'RSV B', ie no filter 
        will be applied to RSV genomes that have been uploaded without an RSV label. This would include all 
        sequences from RefSeq. This should not cause any issues because this filter is only run if we are creating 
        a custom taxonomy for RSV and that means we will not be using any RefSeq or other RSV sequences that 
        have not been loaded into the DB with an RSV label.  But it would make sense to implement a name-based 
        filter that would need a regex to recognize RSV genomes (by full name and acronym).  

        Args:
            categories: list, optional, defaults to ['RSV A','RSV B']
                A list of category labels of sequences to be filtered
                
        Return:
            True on success
            
        Side effects:
            Sets sequences.include values
        """
        # this is here rather than in the method defaults because of this: https://pylint.readthedocs.io/en/latest/user_guide/messages/warning/dangerous-default-value.html
        if not categories:
            categories= ['RSV A','RSV B']
        
        logging.info( f'starting filter to remove incomplete RSV genomes')

        # The RSV genome is a single-stranded, non-segmented molecule that is 15,191–15,226 nucleotides long 
        # https://www.nature.com/articles/s41598-023-40760-y
        # Using a cutoff of 15k
        rsv_min_len = 15000
        n_filtered= 0
        for label in categories:
            sequence_ids_to_remove  = self._db.get_seq_ids_by_category_and_seq_lt(category= label, seq_len_lt= rsv_min_len)
            self._db.mark_as_not_included(sequence_ids_to_remove)
            n_filtered+= len(sequence_ids_to_remove)
        logging.info( f'filtered out {n_filtered} RSV genomes as incomplete (<{rsv_min_len} bases)')
        return True

    def link_all_unlinked_sequences_to_taxonomy_nodes(self):
        """
        For all sequences that do not yet have a tax_id set, this method attempts to set the tax_id using the 
        data in the acc2taxids table, if one was created. The information in the acc2taxids table originates in the 
        NCBI accession-to-taxid file and links GenBank or RefSeq accession IDs to taxonomy IDs.  
        For a sequence to be linked successfully, it must have an NCBI accession ID that was successfully parsed 
        from the FASTA header at the sequences import step and that accession must exist in the NCBI acc2taxid file.  
        Sequences that were obtained from non-NCBI sources may not have an accession assigned and, hence, may not 
        be linked by this method.  
        The tax_id identified in the linkage process is written to the sequences.tax_id field and will therefore be 
        used in the export of sequences to FASTA files.  
        
        Returns:
            True on success
            
        Side-effects:
            Updates sequences.tax_id field for some records
        """
        logging.info('Starting to link currently unlinked sequences to taxonomy')
        n_linked= 0
        seq2taxid_iterator = self._db.all_seq2taxid_iterator(unlinked_only= True)       
        with self._db.bulk_update_buffer(table_name='sequences', id_field='id', update_fields= ['tax_id'], buffer_size= 50000) as seq_update_buffer:
            for row in seq2taxid_iterator:
                n_linked+=1
                n_updated_seqs= seq_update_buffer.add_row(
                    {
                        'id': row['id'],
                        'tax_id': row['tax_id']
                    }
                )
                if n_updated_seqs > 0:
                    logging.info(f'flushed {n_updated_seqs} sequence record updates to DB')
        logging.info(f'Linkage task completed: {n_linked} sequences have been linked to taxonomy nodes' )
        return True

    def create_db_ready_dir( self, path: str, force: bool = True, fasta_file_name:str = 'library.fna' ):
        """
        Create a directory with all files needed to build a new kraken database.
        The FASTA file is named library.fna as per kraken2 convention but can be given a different name. The 
        NCBI taxonomy files names.dmp and nodes.dmp cannot be renamed as kraken2 relies on those names.
        
        Parameters:
            path: str
                path to the directory into which the files will be written.
                If force is not used, the DIR must not exist yet
                
            force: bool, optional, defaults to True
                if true, force overwrite whatever is already in the output dir. If not, will fail
                if directory already exists
                
        Returns:
            True on success
            
        Side effects:
            Writes files to path
            
        """
        library_path = os.path.join( path , 'library' )
        taxonomy_path = os.path.join( path, 'taxonomy' )
        
        fasta_file_path = os.path.join( library_path, fasta_file_name )

        if os.path.exists( path ): 
            if not force:
                raise ValueError(f'directory { path } exists already. Will not write into existing directory')
        else:
            os.makedirs( path )
            
        if not os.path.exists( library_path ):
            os.makedirs( library_path )
        if not os.path.exists( taxonomy_path ):
            os.makedirs( taxonomy_path )
        
        logging.info( f'writing output files to {path}')
        write_fasta(self._db, fasta_file_path)
        write_taxonomy(self._db, taxonomy_path)
        logging.info( 'output files completed')
        return True

    def prune_db(self):
        """
        Remove unnecessary data from the DB at the end of the process.  
        Just delegates to a method of the same name in the Db class. See there for more information.  
        """
        logging.info('starting to prune unnecessary data from the DB')
        self._db.prune_db()
        logging.info('finished DB prune')
        return True