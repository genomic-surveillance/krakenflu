import sqlite3
import os
from collections import defaultdict

"""
This module contains the definitions of classes to handle the sqlite DB for kraken_flu
The schema of the DB is as follows:

  ┌────────────────┐           ┌───────────────────────────────────┐            
  │taxonomy_names  │           │ taxonomy_nodes                    │            
  ├────────────────┤           ├───────────────────────────────────┤            
  │PK id           │    ┌─────►│ PK tax_id                         │▲───┐       
  │FK tax_id       ├────┤      │ FK parent_tax_id                  ├────┘       
  │   name_txt     │    │      │    rank                           │            
  │   unique_name  │    │      │    embl_code                      │            
  │   name_class   │    │      │    division_id                    │            
  └────────────────┘    │      │    inherited_div_flag             │            
                        │      │    genetic_code_id                │            
  ┌──────────────────┐  │      │    inherited_GC_flag              │            
  │sequences         │  │      │    mitochondrial_genetic_code_id  │            
  ├──────────────────┤  │      │    inherited_MGC_flag             │            
  │PK id             │  │      │    GenBank_hidden_flag            │            
  │FK tax_id         ├──┘      │    hidden_subtree_root_flag       │            
  │   fasta_header   │         │    comments                       │            
  │   dna_sequence   │         └───────────────────────────────────┘            
  │   seq_length     │                                                          
  │   ncbi_acc       │                                                          
  │   flu_name       │                                                          
  │   flu_a_h_subtype│                                                          
  │   flu_a_n_subtype│                                                          
  │   include        │                                                          
  └──────────────────┘                                                          

"""

class Db():   
    """
    This class handles all interactions with the SQLite DB, including the DB session.

    Parameters:
        db_path: str, required
            absolute path to the sqlite file to be created
            
        debug: bool, optional, defaults to False
            If True, will emit debug messages from SQLAlchemy
    """
    def __init__( self, db_path: str, debug:bool=False):
        self.db_path = db_path
        if os.path.exists(db_path):
            raise Exception(f"Cannot create database at `{db_path}` - the file exists already.")

        # connect and create a cursor (session)
        self._con = sqlite3.connect(db_path)
        self._con.row_factory = sqlite3.Row # see https://docs.python.org/3.8/library/sqlite3.html?highlight=word#row-objects
        self._cur = self._con.cursor()
        
        # create the DB schema
        self._cur.executescript( self.schema )

    def insert_from_dict(self, table_name:str, data:dict):
        """
        Creates an INSERT statement like:
        INSERT INTO table_name(field1, field2 [...]) VALUES(:field1, :field2 [...])
        which is called with the field values in an execute on the DB cursor

        Args:
            data (dict): dictionary of the data, keys must be valid field names
        """
        field_name_list = [] # for table_name(field list)
        value_name_list = [] # for VALUES(names list)
        values = [] # the values we are inserting in the order of field names
        for key, value in data.items():
            field_name_list.append(key)
            value_name_list.append(':'+key)
            values.append(value)
        stmt = f"INSERT INTO {table_name}({','.join(field_name_list)}) VALUES({','.join(value_name_list)})"
        self._cur.execute(stmt, data)        
        self._con.commit()

    def add_sequence( self, fasta_header:str,  dna_sequence:str, category:str, flu_type:str, ncbi_acc:str, original_taxid:int, is_flu:bool, isolate_name:str, segment_number:int, h_subtype:int, n_subtype:int ):
        """
        Add a sequence record to table "sequences" without a link to a taxon node (which will be provided later).  
        NOTE on tax_id: If the FASTA header contains a kraken:taxid tag, we record this in the original_tax_id field.  
        The tax_id field is only populated later by this tool when we make the final decision about which node we are associating 
        this sequence with. So the original_tax_id is just a hint that may or may not end up being the tax_id we use in the end, but 
        both are recorded.
        """
        self.insert_from_dict('sequences',
            {
                'fasta_header':fasta_header,
                'dna_sequence':dna_sequence,
                'seq_length':len(dna_sequence),
                'ncbi_acc':ncbi_acc,
                'is_flu':int(is_flu),
                'flu_name':isolate_name,
                'flu_type':flu_type,
                'flu_a_h_subtype':h_subtype,
                'flu_a_n_subtype':n_subtype,
                'segment_number':segment_number,
                'include':int(True),
                'category':category,
                'original_tax_id':original_taxid,
            } 
        )
        
    def add_node(self, tax_id:int, parent_tax_id: int, rank: str, embl_code:str, division_id: int, inherited_div_flag:int, genetic_code_id:int,inherited_GC_flag:int, mitochondrial_genetic_code_id:int, inherited_MGC_flag:int, GenBank_hidden_flag:int, hidden_subtree_root_flag:int, comments:str):
        """
        Add a row into the nodes table.  
        """
        self.insert_from_dict('taxonomy_nodes',
            {
                'tax_id': tax_id,
                'parent_tax_id': parent_tax_id,
                'rank': rank,
                'embl_code': embl_code,
                'division_id': division_id,                   
                'inherited_div_flag': inherited_div_flag,            
                'genetic_code_id': genetic_code_id,
                'inherited_GC_flag': inherited_GC_flag,
                'mitochondrial_genetic_code_id': mitochondrial_genetic_code_id,
                'inherited_MGC_flag': inherited_MGC_flag,
                'GenBank_hidden_flag': GenBank_hidden_flag,
                'hidden_subtree_root_flag': hidden_subtree_root_flag,
                'comments': comments
            } 
        )
        
    def get_leaf_node_tax_ids(self):
        """
        Runs a query to fetch all tax_ids of leaf nodes, ie nodes that have no children
        """
        stmt = """
            SELECT tax_id
            FROM taxonomy_nodes AS tn1
            WHERE NOT EXISTS(
                SELECT *
                FROM taxonomy_nodes AS tn2
                WHERE tn2.parent_tax_id = tn1.tax_id
            )
        """
        rows = self._cur.execute(stmt).fetchall()
        return [ r['tax_id'] for r in rows ]

    def add_name(self, tax_id:int, name:str, name_class:str, unique_name:str):
        """
        Add a row into the names table. 
        NOTE: there is no check here to ensure the referenced tax_id actually exists in the taxonomy_nodes 
        parent table. This saves a significant number of SELECTs and the check should not be required because 
        we are basically just re-creating the taxonomy database from NCBI using the flat file export.
        """
        self.insert_from_dict('taxonomy_names',
            {
                'tax_id': tax_id,
                'name': name,
                'name_class': name_class,
                'unique_name': unique_name
            } 
        )
        
    def get_all_tax_ids_paths_root_to_leaf(self):
        """
        Queries the DB to retrieve a list of lists as follows:
            [ 
                [root_tax_id, child_tax_id1, child_tax_id2, ... leaf_tax_id1],
                [root_tax_id, child_tax_id3, child_tax_id4, ... leaf_tax_id2],]
            ]
        Every top level list element is a list of tax_ids that represent a path through the
        taxonomy from root to a leaf. The final datastructure contains one top level list for
        every leaf node in the taxonomy

        Returns:
            list of lists
        """
        leaf_node_tax_ids = self.get_leaf_node_tax_ids()
        data = []
        
        for leaf_node_tax_id in leaf_node_tax_ids:
            data.append( self.get_tax_ids_path_root_to_node(leaf_node_tax_id))
            
        return data

    def get_tax_ids_path_root_to_node(self, starting_tax_id:int):
        """
        For a given node by tax_id, traverse the taxonomy from child to parent and collect 
        tax_ids until we hit the root node.  
        Return as a list in order from root tax_id to the one we started from, ie from parent to 
        child. 
        Example:
        path = db.get_tax_ids_path_root_to(12)
        [1, 4, 12]
        -> means the node with tax_id 12 has a path (of tax_ids) from root (1) via node 4.
        
        Args:
            starting_tax_id: int, required
                The tax_id from where we start to query towards the root.
                This does not have to be a leaf node.

        Returns:
            list: list of tax_ids in order from root to the node we started from
        """
        tax_id = starting_tax_id
        path = []
        while tax_id:
            path.append(tax_id)
            tax_id = self.get_parent_tax_id(tax_id)
        path.reverse()
        return path
        
    
    def get_parent_tax_id(self, tax_id:int):
        """
        From a given node, identified by tax_id, find the parent tax_id.  
        Returns None for the root taxonomy node, which has no tax_id

        Args:
            tax_id: int, required
                The tax_id for the node for which we want to identify the parent

        Returns:
            parent_tax_id: int
        """
        stmt = """
            SELECT tax_id
            FROM taxonomy_nodes
            WHERE tax_id = (
                SELECT parent_tax_id 
                FROM taxonomy_nodes 
                WHERE tax_id = ?
            )
        """
        rows= self._cur.execute(stmt,[tax_id]).fetchall()
        if not rows:
            return None
        elif len(rows) > 1:
            raise ValueError(f"taxonomy node with tax_id {tax_id} has more than one parent, which should not be possible")
        else:
            return rows[0]['tax_id']
        
    def get_flu_segment_data_dict(self):
        """
        Returns a dictionary of all flu names with segment lengths
        Returns:
            Dictionary with the following structure:
                { flu_name: {segment_number: sequence_length} }
        """
        stmt = """
        SELECT 
            flu_name,
            segment_number,
            seq_length
        FROM sequences
        WHERE is_flu = 1
        """
        data = defaultdict(lambda: defaultdict(int))
        rows = self._cur.execute(stmt).fetchall()
        for row in rows:
            data[ row['flu_name'] ][ row['segment_number']] = row['seq_length']
            
        return data
    
    @property        
    def schema(self):
        """
        The database schema as a string. Used for the creation of the DB
        """
        return """
            CREATE TABLE taxonomy_nodes (
                tax_id INTEGER NOT NULL, 
                parent_tax_id INTEGER, 
                rank VARCHAR NOT NULL, 
                embl_code VARCHAR, 
                division_id INTEGER NOT NULL, 
                inherited_div_flag INTEGER NOT NULL, 
                genetic_code_id INTEGER NOT NULL, 
                "inherited_GC_flag" INTEGER NOT NULL, 
                mitochondrial_genetic_code_id INTEGER NOT NULL, 
                "inherited_MGC_flag" INTEGER NOT NULL, 
                "GenBank_hidden_flag" INTEGER NOT NULL, 
                hidden_subtree_root_flag INTEGER NOT NULL, 
                comments VARCHAR, 
                PRIMARY KEY (tax_id), 
                FOREIGN KEY(parent_tax_id) REFERENCES taxonomy_nodes (tax_id)
            );

            CREATE TABLE taxonomy_names (
                id INTEGER NOT NULL, 
                tax_id INTEGER NOT NULL, 
                name VARCHAR NOT NULL, 
                name_class VARCHAR NOT NULL, 
                unique_name VARCHAR, 
                PRIMARY KEY (id), 
                FOREIGN KEY(tax_id) REFERENCES taxonomy_nodes (tax_id)
            );

            CREATE TABLE sequences (
                id INTEGER NOT NULL, 
                tax_id INTEGER, 
                fasta_header VARCHAR NOT NULL, 
                dna_sequence VARCHAR NOT NULL, 
                seq_length INTEGER NOT NULL, 
                segment_number INTEGER, 
                ncbi_acc VARCHAR, 
                flu_name VARCHAR, 
                flu_type VARCHAR, 
                flu_a_h_subtype INTEGER, 
                flu_a_n_subtype INTEGER, 
                include BOOLEAN NOT NULL, 
                is_flu BOOLEAN NOT NULL, 
                category VARCHAR, 
                original_tax_id INTEGER, 
                PRIMARY KEY (id), 
                FOREIGN KEY(tax_id) REFERENCES taxonomy_nodes (tax_id)
            );
        """
