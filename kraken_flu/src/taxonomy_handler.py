import re
import os.path
from cached_property import cached_property
import logging
from kraken_flu.src.utils import parse_flu

logging.basicConfig( format='%(asctime)s %(message)s', level=logging.DEBUG )

class TaxonomyHandler():
    """
    This class handles the NCBI taxonomy, its main task is to restructure the taxonomy for
    influenza isolates.
    
    The flu taxonomy is modified by creating new levels of nodes in order to break down flu references into 
    segments. This creates a structure as shown here for Flu A:                                                                               
                                                                                                            
                                ┌──────────────┐                                                           
                                │ Influenza A  │                                                           
               ┌────────────────┴──────┬───────┴──────────────────────┐                                    
               │                       │                              │                                    
               │                       │                              │                                    
               │                       │                              │                                    
    ┌──────────▼──────────┐ ┌──────────▼──────────┐       ┌───────────▼──────────┐                         
    │Influenza A segment 1│ │Influenza A segment x│ [...] │Influenza A segment 4 │ [...]                   
    └─┬───────────────────┘ └─────────────────────┘  ┌────┴──────────────────┬───┘                         
      │                                              │                       │                             
      │ ┌──────────────────────────┐   ┌─────────────▼──────────┐ ┌──────────▼─────────────┐               
      ├─► isolate segment sequence │   │Influenza A H1 segment 4│ │Influenza A H2 segment 4│ [...]         
      │ └──────────────────────────┘   └┬───────────────────────┘ └──────┬─────────────────┘               
      │                                 │                                │                                 
      │ ┌──────────────────────────┐    │  ┌──────────────────────────┐  │  ┌──────────────────────────┐   
      └─► isolate segment sequence │    ├──► isolate segment sequence │  ├──► isolate segment sequence │   
        └──────────────────────────┘    │  └──────────────────────────┘  │  └──────────────────────────┘   
                                        │                                │                                 
                                        │  ┌──────────────────────────┐  │  ┌──────────────────────────┐   
                                        └──► isolate segment sequence │  └──► isolate segment sequence │   
                                           └──────────────────────────┘     └──────────────────────────┘   
    For Flu B, the same structure is created except that there are no special nodes for the subtypes, ie segments 
    4 and 6 are not getting special treatment,                                                                                                        
    
    Parameters:
        taxonomy_path: str, required
            Path to the taxonomy directory from kraken2-build process (kraken2 default name: taxonomy)
            It needs to contain the names.dmp and nodes.dmp files from NCBI
    """
    
    def __init__( self, taxonomy_path: str ):
        
        if not os.path.isdir( taxonomy_path ):
            raise ValueError(f'path { taxonomy_path} does not exist or is not a directory')
        else:
            self.taxonomy_path = taxonomy_path        
        
        names_file_path = os.path.join(self.taxonomy_path, 'names.dmp')
        if os.path.exists( names_file_path ) and os.path.isfile( names_file_path ):
            self.names_file_path = names_file_path
        else:
            raise ValueError(f'missing file { names_file_path } in taxonomy path {self.taxonomy_path}')

        nodes_file_path = os.path.join(self.taxonomy_path, 'nodes.dmp')
        if os.path.exists( nodes_file_path ) and os.path.isfile( nodes_file_path ):
            self.nodes_file_path = nodes_file_path
        else:
            raise ValueError(f'missing file { nodes_file_path } in taxonomy path {self.taxonomy_path}')
        
        self._nodes = self._read_nodes()
        self._names = self._read_names()
        
        # These will be created if/when influenza/segment taxa are added to the taxonomy
        # Once populated, they return dicts of dicts to map influenza type/subtype/isolate + 
        # segment number to the newly created taxon IDs 
        self._influenza_type_segment_tax_ids = None
        self._influenza_subtype_segment_tax_ids = None
        self._influenza_isolate_segment_tax_ids = None
        
    def _read_nodes(self):
        """
        Populates the "nodes" property by reading the taxonomy file.
        See "nodes" for explanation
        """ 
        logging.info( f'reading NCBI nodes file{ self.nodes_file_path }')
        data = {}
        n = 0
        with open( self.nodes_file_path, 'r' ) as fh:
            for row in fh:
                d = self._read_tax_data_file_row( row )
                tax_id = int(d[0])
                parent_tax_id = int(d[1])
                if tax_id in data:
                    raise ValueError(f'a taxon ID ({tax_id}) was found more than once in file {self.nodes_file_path}, which should not be possible')
                data[ tax_id ] = {
                    'parent_id': parent_tax_id,
                    'data': d[2:]
                }
                n +=1
        logging.info( f'found { n } nodes in { self.nodes_file_path }')
        return data
    
    @property
    def nodes(self):
        """
        An in-memory representation of the contents of the nodes.dmp file as a dictionary with taxon ID as key.
        
        The NCBI nodes.dmp file contains the following data (according to NCBI README):
        tax_id                                  -- node id in GenBank taxonomy database
        parent tax_id                           -- parent node id in GenBank taxonomy database
        rank                                    -- rank of this node (superkingdom, kingdom, ...) 
        embl code                               -- locus-name prefix; not unique
        division id                             -- see division.dmp file
        inherited div flag  (1 or 0)            -- 1 if node inherits division from parent
        genetic code id                         -- see gencode.dmp file
        inherited GC  flag  (1 or 0)            -- 1 if node inherits genetic code from parent
        mitochondrial genetic code id           -- see gencode.dmp file
        inherited MGC flag  (1 or 0)            -- 1 if node inherits mitochondrial gencode from parent
        GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
        hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
        comments                                -- free-text comments and citations

        We only need to manipulate taxon ID and parent taxon ID. The rest of the data is read into a list with the
        same order as shown above, which is also the order in the file.
        
        NOTE: the files is large but not prohibitively large (currently 180MB) so should be ok to just
        slurp into RAM. Makes processing a lot easier.
            
        Returns:
            dict of dicts with the following structure:
                {
                    taxon ID: {
                        'parent_id': parent tax_id ,
                        'data': [ LIST OF ALL REMAINING DATA ITEMS from rank to comments in the above order]
                    }
                }
        """
        return self._nodes
    
    @property
    def names(self):
        """
        An in-memory representation of the contents of the names.dmp file as a dictionary with taxon ID as key. In the case of the names
        file, the taxon ID is not unique across the file. Each taxon ID can have more than one name.
        
        The NCBI names.dmp file contains the following data (according to NCBI README):
        tax_id                                  -- the id of node associated with this name
        name_txt                                -- name itself
        unique name                             -- the unique variant of this name if name not unique
        name class                              -- (synonym, common name, ...)
        
        NOTE: the files is large but not prohibitively large (currently 2300MB) so should be ok to just
        slurp into RAM. Makes processing a lot easier.
            
        Returns:
            dict of list with the following structure:
                {
                    taxon ID: [
                        {
                            'name': name ,
                            'uname': unique name,
                            'nclass': name class (see above)
                        }
                    ]
                }
        """
        return self._names
    
    def _read_names(self):
        """
        Populates the "names" property by reading the taxonomy file.
        See "names" for explanation
        """
        logging.info( f'reading NCBI names file{ self.names_file_path }')

        data = {}
        n = 0
        with open( self.names_file_path, 'r' ) as fh:
            for row in fh:
                d = self._read_tax_data_file_row( row )
                tax_id = int(d[0])
                n += 1
                this_name = {
                    'name': d[1],
                    'uname': d[2],
                    'nclass': d[3]
                }
                if tax_id in data:
                    data[ tax_id ].append( this_name )
                else:
                    data[ tax_id ] = [ this_name ]
        
        logging.info( f'found { n } names in { self.names_file_path }')
        return data
    
    def max_tax_id(self):
        """
        Identifies the maximum (numeric sort) taxon ID in the data.
        This can change as we add more taxa to the taxonomy.
        
        Returns:
            highest tax ID (int)
        """
        max_nodes_tax_ids = sorted(list(self.nodes.keys()))[-1]
        max_names_tax_ids = sorted(list(self.names.keys()))[-1]
        return sorted([max_nodes_tax_ids, max_names_tax_ids])[-1]
    
    def add_taxon( self, tax_id: int, parent_tax_id: int, name: str ):
        """
        Add a new (artificial) taxon to the taxonomy. 
        A new taxon has to have a name, a taxon ID and a parent taxon ID.
        Adding a taxon creates a new entry in "names" and "nodes".
        The parent_tax_id must already exist in "names" and "nodes", so it if is also new,
        it must be inserted first before any children are inserted.
        
        NOTE: The entry in "nodes" is populated with some hardcoded data that makes sense in 
        the context of adding a new taxon to the influenza taxonomy but might need to 
        be looked at if this is to be expanded. 
        
        Parameters:
            tax_id: int, required
                The taxon ID of the new taxon. This needs to be new and cannot exist in the
                taxonomy yet.
                
            parent_tax_id: int, required
                Taxon ID of the parent of this new taxon. Must already exist in the taxonomy.
                
            name: str, required
                Name of the new taxon
                
        Returns:
            True on success
            
        """
        if parent_tax_id not in self.nodes or parent_tax_id not in self.names:
            raise ValueError(f'parent _tax_id {parent_tax_id} does not exist in taxonomy. Create it first, then add children')
        if tax_id in self.nodes or tax_id in self.names:
            raise ValueError(f'tax_id {tax_id} exists already in taxonomy. Must have a new tax_id that does not exist yet.')
        
        self.names[ tax_id ] = [
            {
                'name': name,
                'uname': name,
                'nclass': 'scientific name'
            } 
        ]

        # hardcoded fields: '9' in data field 3 is for viral division. If ever we need to 
        # make this more generic, this is probably the only one that would need to be changed
        self.nodes[ tax_id ] = {
            'parent_id': parent_tax_id,
            'data': [ 'no rank', '', '9', '1', '1', '1', '0', '1', '1', '', '']
        }
        
        logging.info( f'added taxon to taxonomy:  { name }; tax_id: { tax_id }; parent tax_id: { parent_tax_id}')
        
        return True
    
    def write_nodes_file( self, path ):
        """
        Write the nodes data to a new file in the NCBI nodes.dmp format.
        
        Parameters:
            path: str, required
                path to the file that is to be written
                
        Returns:
            True on success
            
        Side-effects:
            Writes to file
            
        """
        logging.info( f'writing nodes file to { path }')

        with open( path, 'w', encoding="utf-8" ) as out_fh:
            for tax_id, node in self.nodes.items():
                out_fh.write(
                    self._format_tax_data_file_output_row(
                        [ tax_id, node['parent_id'] ] + node['data']
                    ) +"\n"
                )
        logging.info( f'finished writing nodes file to { path }')

        return True
        
    def write_names_file( self, path ):
        """
        Write the names data to a new file in the NCBI nodes.dmp format.
        
        Parameters:
            path: str, required
                path to the file that is to be written
                
        Returns:
            True on success
            
        Side-effects:
            Writes to file
            
        """
        logging.info( f'writing names file to { path }')

        with open( path, 'w', encoding="utf-8" ) as out_fh:
            for tax_id, names in self.names.items():
                for name in names:
                    out_fh.write(
                        self._format_tax_data_file_output_row(
                            [ tax_id, name['name'], name['uname'], name['nclass'] ]
                        ) +"\n"
                    )
        logging.info( f'finished writing names file to { path }')

        return True
    
    def create_influenza_type_segment_taxa(self):
        """
        For each of the influenza types for which we want to split genomes into segments (currently 
        A, B and C), one new node is created in the taxonomy for each segment as a child of the 
        "Influenza X" original parent node.
        The parent node ("Influenza A" etc) must exist in the provided taxonomy already. The result of this
        method is a new level of "artificial" taxon nodes like this (for FluA, same for other Flu types):                                        
                                                                                    
                                    ┌──────────────┐                                 
                                    │ Influenza A  │                                 
                   ┌────────────────┴───────┬──────┴───────────────────┐             
        ┌──────────▼──────────┐   ┌─────────▼───────────┐  ┌───────────▼──────────┐  
        │Influenza A segment 1│   │Influenza A segment 2│  │Influenza A segment 3 │  
        └──────-────────────-─┘   └─────────────────────┘  └──────────────────────┘  
                
        NOTE: This is currently restricted to flu A, B and C. If more types should be handled, add the type letter
        to the list "types" and retrieve the respective parent node IDs in teh code below.
        insert the influenza A segment 1, 2, 3 etc nodes. 
        
        The datastructured that is returned can be used later to retrieve the new taxon IDs for the purpose of 
        assigning genome sequences to the new taxon nodes.
        
        Returns:
            sets and returns self.influenza_type_segment_tax_ids
            (datastructure that maps type/segment to the new taxid)
        
        """
        # check if we have already done this and, if so, just return the existing data
        if self.influenza_type_segment_tax_ids is not None:
            return self.influenza_type_segment_tax_ids
        
        # this will be a lookup of parent tax IDs
        flu_tax_ids = {}
        
        # add more types here if needed (e.g. 'B') - would need the parent node IDs as well (above)
        types = ['A','B','C','D']
        data = {}
        new_tax_id_i = self.max_tax_id()
        for type in types:
            flu_tax_ids[type], _ = self._tax_id_and_parent_id_by_name( f'Influenza {type} virus')
            if flu_tax_ids[type] is None:
                logging.info( f'could not find Influenza {type} in names file: not creating taxonomy for this influenza type')
                continue
            for seg_num in range(1,9):
                new_tax_id_i += 1
                self.add_taxon( tax_id= new_tax_id_i, parent_tax_id= flu_tax_ids[type], name=f'Influenza {type} segment {seg_num}' )
                try:
                    data[type][seg_num] = new_tax_id_i
                except KeyError:
                    data[type] = { seg_num: new_tax_id_i}
            
        self._influenza_type_segment_tax_ids = data
        return data
    
    @property
    def influenza_type_segment_tax_ids(self):
        """
        Datastructure of taxon IDs for influenza type + segment. Populated when running 
        self.create_influenza_type_segment_taxa
        
        Returns:
            Dict of dicts with the following structure:
            
                {
                    TYPE: {
                        SEGMENT_NUMBER: tax_id
                    }
                }
                
        """
        return self._influenza_type_segment_tax_ids
            
    
    def create_influenza_subtype_segment_taxa(self):
        """
        Similar to "create_influenza_type_segment_taxa". Inserts new taxa into the taxonomy for 
        Influenza subtype/segment number such as "Influenza A H1 segment 4", "Influenza A N3 segment 6"
        
        In this current implementation, we only insert this level of the taxonomy for segments 4 and 6 of 
        influenza A, because only Influenza A subtyping is built on these segments (Ha and Na genes).
        For Influenza B, no additional level is inserted here and, consequently, all Influenza B segment sequences 
        are children of the taxa created in create_influenza_type_segment_taxa. 
        For segment 4, we create one taxon for each H subtype (segment 4 encodes the HA gene). For
        segment 6, one taxon is added for each N subtype (6 encodes NA gene).
        
        NOTE: in the current implementation, the 'segment' level in the resulting data structure is
        redundant because every H subtype only has segment 4 entries and every N subtype only has segment 6.
        Currently keeping that level in the data structure simply because it might be used in the future.
        
        There are 18 known H subtypes and 11 known N subtypes. We create all 18 "Hx segment 4" nodes as well
        as all 11 "Nx segment 6" nodes.
        
        The new nodes are child nodes to the type/segment nodes. If create_influenza_type_segment_taxa 
        has not been run, it is run here.
        
        Returns:
            sets and returns self.influenza_subtype_segment_tax_ids
        
        """
        # check if we have already done this and, if so, just return the existing data
        if self.influenza_subtype_segment_tax_ids is not None:
            return self.influenza_subtype_segment_tax_ids
        
        if not self.influenza_type_segment_tax_ids:
            self.create_influenza_type_segment_taxa()
        
        data={}
            
        # create the Hx segment 4 nodes as children of Influenza A segment 4 node
        parent_tax_id = self.influenza_type_segment_tax_ids['A'][4]
        if not parent_tax_id:
            raise ValueError('cannot retrieve taxon ID for Influenza A segment 4')
        
        new_tax_id_i = self.max_tax_id()
        for subtype_num in range(1,19):
            new_tax_id_i+=1
            subtype = 'H'+str(subtype_num)
            self.add_taxon( tax_id=new_tax_id_i, parent_tax_id=parent_tax_id, name=f'Influenza A {subtype} segment 4')
            try:
                data[subtype][4] = new_tax_id_i
            except KeyError:
                data[subtype] = { 4: new_tax_id_i}
            
        # create the Nx segment 4 nodes as children of Influenza A segment 6 node
        parent_tax_id = self.influenza_type_segment_tax_ids['A'][6]
        if not parent_tax_id:
            raise ValueError('cannot retrieve taxon ID for Influenza A segment 6')
        
        for subtype_num in range(1,12):
            new_tax_id_i+=1
            subtype = 'N'+str(subtype_num)
            self.add_taxon( tax_id=new_tax_id_i, parent_tax_id=parent_tax_id, name=f'Influenza A {subtype} segment 6')
            try:
                data[subtype][6] = new_tax_id_i
            except KeyError:
                data[subtype] = { 6: new_tax_id_i}

        self._influenza_subtype_segment_tax_ids = data
        return data
    
    @property
    def influenza_subtype_segment_tax_ids(self):
        """
        Datastructure of taxon IDs for influenza subtype + segment. Populated when running 
        self.create_influenza_subtype_segment_taxa       
        
        Returns:
            Dict of dicts with the following structure:
            
                {
                    SUBTYPE: {
                        SEGMENT_NUMBER: tax_id
                    }
                }
        """
        return self._influenza_subtype_segment_tax_ids
    

    def create_influenza_isolate_segment_taxa(self):
        """
        Like self.create_influenza_type_segment_taxa and self.create_influenza_subtype_segment_taxa,
        this method creates new nodes in the taxonomy. The nodes created in this method represent the actual 
        reference genomes, ie the isolate segment sequences. These are added to the new levels of the 
        taxonomy created by self.create_influenza_type_segment_taxa and self.create_influenza_subtype_segment_taxa.
        
        For each influenza isolate, 8 new nodes are created in the taxonomy, one for each of the possible 
        segments of the isolate. This class does not look at the available sequence data, thus is does 
        not matter whether or not genome sequence exists for this isolate for all 8 segments. The same 
        applies to all nodes of the taxonomy. They exist whether or not we have a sequence. This follows the general 
        design of the NCBI taxonomy files, which can be linked to but do not depend on sequence records.
        
        Influenza A and B/C are treated differently because of the extra level in the taxonomy for segments 4 
        and 6 of influenza A, which we don't for influenza B or C.
        
        Returns:
            sets and returns self.influenza_isolate_segment_tax_ids
        """
        # check if we have already done this and, if so, just return the existing data
        if self.influenza_isolate_segment_tax_ids is not None:
            return self.influenza_isolate_segment_tax_ids
        
        # we need type/segment and subtype/segment nodes already done, make sure this is the case
        # or trigger that cascade now
        if not self.influenza_subtype_segment_tax_ids:
            self.create_influenza_subtype_segment_taxa()
        
        data = {}
        new_tax_id_i = self.max_tax_id() + 1
        # identify the isolates of interest
        # which is currently defined by nodes with a name of a flu A isolate
        # take a snapshot of keys before we start inserting new ones otherwise we would be
        # iterating over a dict with changing size as we are adding new nodes
        existing_tax_ids = [x for x in self.names.keys()]
        for tax_id in existing_tax_ids:
            for name_record in self.names[ tax_id ]:
                name = name_record['name']
                nclass = name_record['nclass'] # scientific name, equivalent name etc (name class)
                
                is_flu, flu_type, isolate_name, h_subtype, n_subtype, _ = parse_flu( name )

                if nclass == 'scientific name' and isolate_name is not None and flu_type is not None:
                    # this is a flu isolate, generate 8 nodes in the taxonomy,
                    # one for each segment. 
                    # Segments 4 and 6 need special treatment because, in flu A, we attach these as 
                    # children of the Influenza A Hx segment 4 and Influenza A Nx segment 6 parent nodes, 
                    # which don't exist for fluB
                    seg4_parent_id = seg6_parent_id = new_tax_name = parent_tax_id = None
                    if h_subtype:
                        try:
                            seg4_parent_id = self.influenza_subtype_segment_tax_ids[h_subtype][4]
                        except KeyError:
                            raise ValueError(f'isolate has H subtype ({h_subtype}) but no parent tax id could be found for segment 4 in this subtype')
                    if n_subtype:
                        try:
                            seg6_parent_id = self.influenza_subtype_segment_tax_ids[n_subtype][6]
                        except KeyError:
                            raise ValueError(f'isolate has N subtype ({n_subtype}) but no parent tax id could be found for segment 6 in this subtype')
                    
                    # create the new taxa
                    for seg_num in range(1,9):
                        new_tax_name = isolate_name + ' segment ' + str(seg_num)
                        if seg_num == 4 and seg4_parent_id:
                            parent_tax_id = seg4_parent_id
                        elif seg_num == 6 and seg6_parent_id:
                            parent_tax_id = seg6_parent_id
                        else:
                            try:
                                parent_tax_id = self.influenza_type_segment_tax_ids[flu_type][seg_num]
                            except KeyError:
                                raise ValueError(f'no parent tax id for flu {flu_type} segment {seg_num}')
                        self.add_taxon(tax_id= new_tax_id_i, parent_tax_id= parent_tax_id, name= new_tax_name)
                        try:
                            data[isolate_name][seg_num] = new_tax_id_i
                        except KeyError:
                            data[isolate_name] = { seg_num: new_tax_id_i }
                        
                        new_tax_id_i += 1

        self._influenza_isolate_segment_tax_ids = data
        return data
    
    @property
    def influenza_isolate_segment_tax_ids(self):
        """
        Datastructure of taxon IDs for influenza subtype + segment. Populated when running 
        self.create_influenza_isolate_segment_taxa

        Dict of dicts with the following structure:
            
                {
                    NAME: {
                        SEGMENT_NUMBER: tax_id
                    }
                }
                
                where NAME is an isolate name such as 'A/California/07/2009(H1N1)'
        """
        return self._influenza_isolate_segment_tax_ids
    
    
    def _read_tax_data_file_row( self, row ):
        """
        Parses one row of data from names and nodes dmp file and returns as list
        Removes the trailing \t| from the last column
        """
        data = row.rstrip().split("\t|\t")
        data[-1] = data[-1].rstrip("\t|")
        return data
    
    def _format_tax_data_file_output_row( self, data:list ):
        """
        Formats a list of data into the NCBI dmp file format. Integers are converted to string.
        
        Parameters:
            data: list, required
                list of data
                
        Returns:
            string - formatted data
        """
        return "\t|\t".join([ str(x) for x in data ]) + "\t|"
    
    def _tax_id_and_parent_id_by_name(self, name:str, nclass:str='scientific name'):
        """
        Finds a node in the "names" by name and name-class. Returns taxon and parent IDs
        
        Parameters:
            name: str, required
                name to search for
                
            nclass: str, optional
                name class, defaults to 'scientific name'
                
        Returns:
            taxon ID, parent taxon ID 
        """
        if name=='':
            return None
        
        found = False
        tax_id_final = None
        parent_id_final = None
        for tax_id, name_records in self.names.items():
            for name_record in name_records:
                if name_record['name'] == name and name_record['nclass'] == nclass:
                    parent_id = self.nodes[tax_id]['parent_id']
                    if not parent_id:
                        raise ValueError(f'could not find a taxon ID in nodes that should be there:{tax_id}')
                    else:
                        if found:
                            raise ValueError(f'found more than one match to name:"{name }" name-class:"{nclass}')
                        else:
                            found = True
                            tax_id_final = tax_id
                            parent_id_final = parent_id       
        return tax_id_final, parent_id_final

