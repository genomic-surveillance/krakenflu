import re
import os.path
from cached_property import cached_property
import logging

logging.basicConfig( format='%(asctime)s %(message)s', level=logging.DEBUG )

class TaxonomyHandler():
    """
    This class handles the NCBI taxonomy, its main task is to restructure the taxonomy for
    influenza isolates.
    
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
        
        
    @cached_property
    def nodes(self):
        """
        Reads the nodes.dmp file into memory as a dictionary with taxon ID as key.
        
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
        data = {}
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
            
        return data
    
    @cached_property
    def names(self):
        """
        Reads the names.dmp file into memory as a dictionary with taxon ID as key. In the case of the names
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
        data = {}
        with open( self.names_file_path, 'r' ) as fh:
            for row in fh:
                d = self._read_tax_data_file_row( row )
                tax_id = int(d[0])
                this_name = {
                    'name': d[1],
                    'uname': d[2],
                    'nclass': d[3]
                }
                if tax_id in data:
                    data[ tax_id ].append( this_name )
                else:
                    data[ tax_id ] = [ this_name ]
            
        return data
    
    def max_tax_id(self):
        """
        Identifies the maximum (numeric sort) taxon ID in the data.
        
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
        
        return True
    
    def _read_tax_data_file_row( self, row ):
        """
        Parses one row of data from names and nodes dmp file and returns as list
        Removes the trailing \t| from the last column
        """
        data = row.rstrip().split("\t|\t")
        data[-1] = data[-1].rstrip("\t|")
        return data