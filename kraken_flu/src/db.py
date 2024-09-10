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

         # Enable foreign key support
        self._con.execute("PRAGMA foreign_keys = ON;")
        
        if debug:
            self._con.set_trace_callback(print)
        
        self._con.row_factory = sqlite3.Row # see https://docs.python.org/3.8/library/sqlite3.html?highlight=word#row-objects
        self._cur = self._con.cursor()
        
        # create the DB schema
        self._cur.executescript( self.schema )

    def insert_from_dict(self, table_name:str, data:dict):
        """
        Creates an INSERT statement like:
        INSERT INTO table_name(field1, field2 [...]) VALUES(:field1, :field2 [...])
        which is called with the field values in an execute on the DB cursor
        
        TODO: this method and bulk_insert are very similar and should be merged into one method
        that can handle both, a single dict or the bulk data list of lists used for the bulk_insert method.  
        This could be done either by checking the type of input data we receive or by explicitly 
        setting a mode argument. The only real difference between the two is how to construct the 
        field values list for the sql. We can use executemany for both use cases.  

        Args:
            data (dict): dictionary of the data, keys must be valid field names
            
        Returns:
            True on success
            
        Side effects:
            Insert data into table
            
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
        return True
        
    def bulk_insert(self, table_name:str, field_names:list, field_data:list):
        """
        Performs a bulk insert of several rows into a database table.  
        
        Args:
            table_name: str, required
                Name of the table to insert into
                
            field_names: list, required
                List of field names to insert into
                
            field_data: list, required
                List of lists of data for the table fields/columns. Order of data in lists must match 
                the order of field provided in argument field_names
                
        Returns:
            True on success
            
        Side effects:
            Insert data into table
            
        """
        if not isinstance(field_data,list) or not isinstance(field_data[0],list):
            raise ValueError("field_data must be a list of lists")
        
        stmt = f"INSERT INTO {table_name}({','.join(field_names)}) VALUES({','.join(['?']*len(field_names))})"
        self._cur.executemany(stmt, field_data)        
        self._con.commit()
        return True
        
    def bulk_insert_buffer(self, table_name:str, buffer_size= 5000):
        """
        Returns a context manager BulkInsertBuffer object to manage buffered bulk inserts.  
        See class definition for BulkInsertBuffer in this module for details.
        
        Args:
            table_name: str, required
                The name of the table we are inserting into
                
            buffer_size: int, optional, defaults to 5000
                The number of rows of data to hold in RAM before flushing to the DB
            
        """
        return BulkInsertBuffer( db=self, table_name= table_name, buffer_size= buffer_size )
        
    def bulk_update(self, table_name:str, update_fields:list, id_field:str, field_data:list, id_field_values:list):
        """
        Performs a bulk update of several rows in a database table. 
        Requires a list of list of update values and a list of ID field values to identify the rows 
        to be updated.  
        To illustrate how the params are used:
        The statement for executemany is constructed as follows:
        
            UPDATE {table_name} SET {update_fields[0]} = ?, {update_fields[1]} = ? [...] WHERE {id_field} = ?
        
        The values for parameter substitution are taken from the two lists field_data and id_field_values.  
        
        Args:
            table_name: str, required
                Name of the table to insert into
                
            update_fields: list, required
                List of field names to update. Must be in the same order as the order of 
                field values in the field_data list of lists.  
                
            id_field: str, required
                Name of the "id_field", which is used in the WHERE clause of the UPDATE statement to 
                identify each row to be updated. Usually the PK but can be any field that can be used in a 
                UPDATER statement WHERE clause.  
                
            field_data: list, required
                List of lists of data for the table fields/columns. Order of data in lists must match 
                the order of field provided in argument field_names
                
            id_field_values: list, required
                List of values for the "id_field", which identifies the rows to be updated
                
        Returns:
            True on success
            
        Side effects:
            Insert data into table
            
        """
        if not isinstance(field_data,list) or not isinstance(field_data[0],list):
            raise ValueError("field_data must be a list of lists")
        
        # combine field update values and the WHERE clause ID field value into a single list
        # with the id field value always the last item in the list
        data = [ [*field_values,id_value] for field_values,id_value in zip(field_data, id_field_values)]

        # build the SQL of the form
        # UPDATE table_name SET field1 = ?, field2 = ? [...] WHERE id_field = ?
        stmt = ' '.join( [
            'UPDATE',
            table_name,
            'SET',
            ', '.join([ f"{x} = ?" for x in update_fields ]),
            'WHERE',
            id_field,
            '= ?'
        ])

        self._cur.executemany(stmt, data)        
        self._con.commit()
        return True

    def bulk_update_buffer(self, table_name:str, id_field:str, update_fields:list, buffer_size= 5000):
        """
        Returns a context manager BulkUpdateBuffer object to manage buffered bulk updates.  
        See class definition for BulkUpdateBuffer in this module for details.
        
        Args:
        
        table_name: str, required
            The name of the table we are inserting into
            
        id_field: str, required
            Name of the field that is used to identify the row that is to be updated.  
            Usually the primary key of the table but doesn't have to be.
            Each row of data has to provide a value for this field.  
            
        update_fields: list(str), required
            List of the field names to update. All rows added to the buffer must provide 
            values for all fields in this list, ie every update must be for the same set of 
            field names. Do not include the name of the ID field in this list, it needs to be 
            se in param "id_field" instead.  
            
        buffer_size: int, optional, defaults to 5000
            The number of rows of data to hold in RAM before flushing to the DB
            
        """
        return BulkUpdateBuffer( db=self, table_name= table_name, id_field= id_field, update_fields= update_fields, buffer_size= buffer_size )

    def add_seq2taxid(self, ncbi_acc:str, acc_version:int, tax_id:int, gi:int,):
        """
        Add a sequence record into the "seq2taxid" table
        """
        self.insert_from_dict('seq2taxid',
            {
                "ncbi_acc":ncbi_acc,
                "acc_version":acc_version,
                "tax_id":tax_id,
                "gi":gi
            }
        )

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
        
    def add_taxon(self, tax_id:int, parent_tax_id:int, name:str):
        """
        TODO: not used in assign_flu_taxonomy_nodes now - is this still needed or can it be deleted?
              Might still be useful in the future for small scale inserts but too slow for bulk operations.  
        Adds a new taxon to the DB, which consists of a new record in taxonomy_nodes and a linked 
        record in taxonomy_names.  
        This is meant to be used for the creation of new "artificial" taxa. Some field data needs to be 
        provided for a node record that is not relevant for these new taxon and is just hardcode here.  
        The name_class is hardcoded to be "scientific name".  
        A parent_tax_id is required so that the new node has a parent in the taxonomy.  

        Args:
            tax_id: int, required
                The tax_id for the new node
        
            parent_tax_id: int, required
                The tax_id of the parent node of this new taxon.  
                
            name: str, required
                The scientific name for the new taxon
        """
        self.add_node(
            tax_id= tax_id,
            parent_tax_id= parent_tax_id,
            rank='no rank',
            embl_code=None,
            division_id=9,
            inherited_div_flag=1,
            genetic_code_id=1,
            inherited_GC_flag=1,
            mitochondrial_genetic_code_id=0,
            inherited_MGC_flag=1,
            GenBank_hidden_flag=0,
            hidden_subtree_root_flag=0,
            comments='kraken_flu added node'
        )
        self.add_name(
            tax_id= tax_id,
            name= name,
            name_class= 'scientific name',
            unique_name= name
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
        Returns None for the root taxonomy node, which has no parent.  
        NOTE: in the conventional adjacency list pattern, the root node should 
        have NULL in the parent_tax_id but the NCBI taxonomy root node has
        tax_id=1 AND parent_tax_id=1
        We generalise this here and assume it is the root node when the node
        refers to itself as parent, ie tax_id==parent_tax_id

        Args:
            tax_id: int, required
                The tax_id for the node for which we want to identify the parent

        Returns:
            parent_tax_id: int
        """
        stmt = """
                SELECT parent_tax_id 
                FROM taxonomy_nodes 
                WHERE tax_id = ?
        """
        rows= self._cur.execute(stmt,[tax_id]).fetchall()
        if not rows:
            return None
        elif len(rows) > 1:
            raise ValueError(f"taxonomy node with tax_id {tax_id} has more than one parent, which should not be possible")
        else:
            parent_tax_id = rows[0]['parent_tax_id']
            if not parent_tax_id or parent_tax_id == tax_id:
                # referring to itself as parent is the hallmark of the root node in NCBI taxonomy
                return None
            else:
                return parent_tax_id
        
    def get_flu_name_segment_data_dict(self):
        """
        Returns a dictionary of all flu names with segment lengths for all records in sequences that
        are marked as flu and have a flu name (there are also records that are marked as flu and have 
        no flu name because it could not be parsed - we can't use those because we can't group segments)
        The returned data structure groups the data by flu isolate name.
        
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
        WHERE is_flu = 1 AND flu_name IS NOT NULL AND segment_number IS NOT NULL
        """
        data = defaultdict(lambda: defaultdict(int))
        rows = self._cur.execute(stmt).fetchall()
        for row in rows:
            data[ row['flu_name'] ][ row['segment_number']] = row['seq_length']
            
        return data
    
    def retrieve_sequence_ids_by_flu_name(self, name:str):
        """
        Query sequences by flu_name and return a list of the sequences.id

        Args:
            name: str, required
                flu_name field value to search with
                
        Returns:
            List of sequences.id for records with matching flu_name value
        """
        stmt = """
            SELECT id
            FROM sequences
            WHERE flu_name = ?
        """
        return [x['id'] for x in self._cur.execute(stmt,[name]).fetchall()]
    
    def retrieve_tax_id_by_node_scientific_name(self,name:str, one:bool=True):
        """
        Retrieves the tax_id for a node identified by its scientific name.  

        Args:
            name: str, required
                Scientific name to search for.  
                
            one: bool, optional, defaults to True
                If True, returns a single tax_id and throws exception if more than 
                one is found
                
        Returns:
            tax_id (int) or list of tax_ids if one==False
            None if no rows found
        """
        stmt = """
        SELECT taxonomy_nodes.tax_id
        FROM taxonomy_nodes
        INNER JOIN taxonomy_names ON(taxonomy_names.tax_id = taxonomy_nodes.tax_id) 
        WHERE name_class = "scientific name" AND name = ? 
        """
        rows = self._cur.execute(stmt,[name]).fetchall()
        if not rows:
            return None
        tax_ids = [x['tax_id'] for x in rows]
        if one:
            if len(tax_ids)>1:
                raise ValueError(f"more than one record found with scientific name '{name}")
            else:
                return tax_ids[0]
        else:
            return tax_ids
    
    def retrieve_unnamed_unsegmented_flu(self):
        """
        Retrieve records from sequences where flu_name is empty or flu segment is not assign but is_flu is True, ie the record
        name identifies it as flu but it doesn't have a proper flu isolate name.
        """
        stmt = """
        SELECT id
        FROM sequences
        WHERE is_flu = 1 AND (flu_name IS NULL or segment_number IS NULL)
        """
        return [ x['id'] for x in self._cur.execute(stmt)]
    
    def retrieve_all_flu_sequences(self, skip_excluded:bool=False):
        """
        Retrieve all sequences marked as flu.
        
        Args:
            skip_excluded: bool, optional, defaults to False
                adds a filter "include=1" to the query
        
        Returns:
            list of rows
        """
        stmt="""
        SELECT 
            id,
            flu_name,
            segment_number,
            flu_type, 
            flu_a_h_subtype, 
            flu_a_n_subtype
        FROM sequences
        WHERE is_flu = 1
        """
        if skip_excluded:
            stmt+= ' AND include = 1'
        return self._cur.execute(stmt).fetchall()
    
    def mark_as_not_included(self, ids:list):
        """
        Marks sequences as "not included" for final output by setting the include value to 0

        Args:
            ids: list, required
                List of sequences.id for records we want to exclude from final output

        Returns:
            True on success
        
        Side effects:
            Sets value of sequences.include
        """
        stmt="""
            UPDATE sequences
            SET include = 0
            WHERE id IN (%s)  
        """
        ids_str = ','.join(map(str,ids))
        self._cur.execute(stmt % ids_str)
        self._con.commit()
        
    def set_tax_id_for_sequence(self, id:int, tax_id:int):
        """
        Set the tax_id for a sequence record identified by its sequence.id  
        
        Args:
            id: int, required
                id of the sequences record to be updated
                
            tax_id: int, required
                new tax_id for the sequences record
                
        Returns:
            True on success
            
        Side effects:
            Sets sequences.tax_id
        """
        stmt="""
            UPDATE sequences
            SET tax_id = ?
            WHERE id= ?  
        """
        self._cur.execute(stmt, [tax_id, id])
        self._con.commit() 
        return True
        
    def set_tax_id_and_mod_fasta_header_for_sequence(self, id:int, tax_id:int, mod_fasta_header:str):
        """
        Set the tax_id for a sequence record identified by its sequence.id  
        TODO: check if this is still needed (should be superseded by bulk_update)
        Args:
            id: int, required
                id of the sequences record to be updated
                
            tax_id: int, required
                new tax_id for the sequences record
                
            mod_fasta_header: str, required
                A modified FASTA header for this record, which will be used when the 
                data is exported to a FASTA file
                
        Returns:
            True on success
            
        Side effects:
            Sets sequences.tax_id
        """
        stmt="""
            UPDATE sequences
            SET tax_id = ?, mod_fasta_header = ?
            WHERE id= ?  
        """
        self._cur.execute(stmt, [tax_id, mod_fasta_header, id])
        self._con.commit() 
        return True
        
    def max_tax_id(self):
        """
        Returns the maximum tax_id from the nodes table, which is used to assign safe new taxon ids for the 
        nodes we are going to add.

        Returns:
            int: maximum tax_id in the nodes table
        """
        stmt="""
        SELECT MAX(tax_id) AS max_tax_id
        FROM taxonomy_nodes
        """
        row = self._cur.execute(stmt).fetchone()
        return row['max_tax_id']
    
    def all_sequences_iterator(self, included_only:bool=True):
        """
        Creates an iterator over all records in the sequences table. By default, filtering 
        for records that have the "include" flag set.  
        
        Args:
            included_only: bool, optional, defaults to True.
                If True, retrieves only records with include=1

        Returns:
            iterator over sequences records.  
            Can be used like this:
                it = db.all_sequences_iterator()
                for row in it:
                    do something with row
        """
        stmt="""
        SELECT
            id, 
            tax_id, 
            fasta_header,
            mod_fasta_header,
            dna_sequence, 
            seq_length, 
            segment_number, 
            ncbi_acc, 
            flu_name, 
            flu_type, 
            flu_a_h_subtype, 
            flu_a_n_subtype, 
            include, 
            is_flu, 
            category, 
            original_tax_id        
        FROM sequences
        """
        if included_only:
            stmt+=' WHERE include = 1'
            
        return self._iterator(stmt)
            
    def all_taxonomy_names_iterator(self):
        """
        Creates an iterator over all records in the names table. 

        Returns:
            iterator over taxonomy_names records.  
        """
        stmt="""
        SELECT
            id, 
            tax_id, 
            name, 
            name_class, 
            unique_name    
        FROM taxonomy_names
        """
        return self._iterator(stmt)
    
    def all_taxonomy_nodes_iterator(self):
        """
        Creates an iterator over all records in the nodes table. 

        Returns:
            iterator over taxonomy_nodes records.  
        """
        stmt="""
        SELECT
            tax_id, 
            parent_tax_id, 
            rank, 
            embl_code, 
            division_id, 
            inherited_div_flag, 
            genetic_code_id, 
            "inherited_GC_flag", 
            mitochondrial_genetic_code_id, 
            "inherited_MGC_flag", 
            "GenBank_hidden_flag", 
            hidden_subtree_root_flag, 
            comments 
        FROM taxonomy_nodes
        """
        return self._iterator(stmt)

    
    def _iterator(self, stmt):
        for row in self._con.execute(stmt):
            yield row
    
    def get_field_names(self, table_name:str):
        """
        Returns a list of the column/field names of the table

        Args:
            table_name: str, required

        Returns:
            list of field names
        """
        data= self._cur.execute(f"SELECT * FROM {table_name} LIMIT 1")
        return [c[0] for c in data.description]
    
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
            CREATE INDEX idx_parent_tax_id
                ON taxonomy_nodes (parent_tax_id);

            CREATE TABLE taxonomy_names (
                id INTEGER NOT NULL, 
                tax_id INTEGER NOT NULL, 
                name VARCHAR NOT NULL, 
                name_class VARCHAR NOT NULL, 
                unique_name VARCHAR, 
                PRIMARY KEY (id), 
                FOREIGN KEY(tax_id) REFERENCES taxonomy_nodes (tax_id)
            );
            CREATE INDEX idx_tax_name 
                ON taxonomy_names (name);
            CREATE INDEX idx_tax_id 
                ON taxonomy_names (tax_id);

            CREATE TABLE sequences (
                id INTEGER NOT NULL, 
                tax_id INTEGER, 
                fasta_header VARCHAR NOT NULL,
                mod_fasta_header VARCHAR NULL,
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

            CREATE TABLE seq2taxid (
                id INTEGER NOT NULL,
                ncbi_acc VARCHAR,
                acc_version INTEGER,
                tax_id INTEGER,
                gi INTEGER,
                PRIMARY KEY (id),
                FOREIGN KEY(ncbi_acc) REFERENCES sequences (ncbi_acc),
                FOREIGN KEY(tax_id) REFERENCES taxonomy_nodes (tax_id)
            );
            CREATE INDEX idx_seq_flu_name 
                ON sequences (flu_name);
            CREATE INDEX idx_seq_tax_id 
                ON sequences (tax_id);
        """

class BulkInsertBuffer():
    """
    This class provides a buffer for INSERT statements into a single table, with the 
    purpose of increasing performance at the database build step where we have to load 
    GB of data into the database from NCBI taxonomy files and FASTA files. This would take 
    too long with one INSERT transaction per row of incoming data.   
    The buffer can be filled with data row by row and automatically triggers a write 
    operation when it is full. The size of the buffer can be set to a desired number of rows.  
    
    Use with context manager like so:
    
        >>> row_data = [{'name':'some name'},{'name':'some other name'}]
        >>> with BulkInsertBuffer('taxonomy_names') as b:
        >>>     for row in row_data:
        >>>         b.add_row(row)
            
    At this point, all row data will have been committed to the DB, nothing else to do
        
    Args:
        table_name: str, required
            The name of the table we are inserting into
            
        buffer_size: int, optional, defaults to 5000
            The number of rows of data to hold in RAM before flushing to the DB
        
    """
    def __init__(self, table_name:str, db:Db, buffer_size= 5000):
        self._db = db
        self.table_name = table_name
        self.field_names = db.get_field_names(table_name)
        self.buffer_size = buffer_size
        self.field_data = [] # a list of lists of row data populated when calling add_row
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, exc_tb):
        if self.field_data:
            self._write_buffer()
            
    def add_row(self, data_dict:dict):
        """
        Adds a row of data to the buffer. When the buffer is full, triggers a DB INSERT. 
        
        Args:
            dict of field names to field values
            
        Return:
            number of flushed rows to DB, 0 when we are just buffering >0 when commit triggered

        """
        row_data = []
        n=0
        for field_name in self.field_names:
            if field_name in data_dict:
                row_data.append( data_dict[field_name])
                n+=1
            else:
                # no data provided for this field, set it to NULL
                row_data.append(None)

        if n < len(data_dict):
            # if we are here, we have unused items in the kwargs dict, ie values for fields that don't exist
            raise ValueError("field names provided to 'add_row' that do not exist in the database table")
        
        self.field_data.append( row_data )
        n_inserted = 0
        if len( self.field_data ) >= self.buffer_size:
            n_inserted = self._write_buffer()

        return n_inserted
        
    def _write_buffer(self):
        """
        Commits the buffer to the database and empties it out. Called by add_row whenever we hit the 
        max buffer size of rows and by __exit__ to empty out the remaining buffer before we lose the 
        context.  
        Returns number of rows inserted
        """
        n = len(self.field_data)
        self._db.bulk_insert(self.table_name, self.field_names, self.field_data)
        self.field_data = [] # reset the buffer
        return n
    
class BulkUpdateBuffer():
    """
    Same as BulkInsertBuffer but for handling bulk UPDATE instead of INSERT operations.  
    TODO: is it worth creating a common base class for this class and BulkInsertBuffer? 
    
    Use with context manager like so:
    
        >>> row_data = [{'tax_id': 123, 'name':'a new name'} ,{'tax_id': 345, 'name':'new name 2'}]
        >>> with BulkUpdateBuffer(table_name='taxonomy_names', id_field='tax_id') as b:
        >>>     for row in row_data:
        >>>         b.add_row(row)
            
    At this point, all row data will have been committed to the DB, nothing else to do
        
    Args:
        table_name: str, required
            The name of the table we are inserting into
            
        id_field: str, required
            Name of the field that is used to identify the row that is to be updated.  
            Usually the primary key of the table but doesn't have to be.
            Each row of data has to provide a value for this field.  
            
        update_fields: list(str), required
            List of the field names to update. All rows added to the buffer must provide 
            values for all fields in this list, ie every update must be for the same set of 
            field names. Do not include the name of the ID field in this list, it needs to be 
            se in param "id_field" instead.  
            
        buffer_size: int, optional, defaults to 5000
            The number of rows of data to hold in RAM before flushing to the DB
        
    """
    def __init__(self, table_name:str, db:Db, id_field:str, update_fields:list, buffer_size= 5000):
        self._db = db
        self.table_name = table_name
        self.id_field = id_field
        self.update_fields = update_fields
        self.buffer_size = buffer_size
        self.field_data = [] # a list of lists of row data populated when calling add_row
        self.id_field_values = [] # list of the ID field values, used to identify the rows to update
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, exc_tb):
        if self.field_data:
            self._write_buffer()
            
    def add_row(self, data_dict:dict):
        """
        Adds a row of data to the buffer. When the buffer is full, triggers a DB UPDATE.
        Keys are field names to be updates EXCEPT the id_field, which is used to identify the
        row to be updated. A value must always be provided for the id_field
        
        Args:
            dict of field names to field values
            
        Return:
            number of flushed rows to DB, 0 when we are just buffering >0 when commit triggered

        """
        row_data = []
        n=0
        try:
            id_value = data_dict[ self.id_field ]
            n+=1
        except KeyError:
            raise ValueError(f"missing value for id_field '{self.id_field}'")
        
        for field_name in self.update_fields:
            if field_name in data_dict:
                row_data.append( data_dict[field_name])
                n+=1
            else:
                raise ValueError(f"no data provided for field {field_name}")

        if n < len(data_dict):
            # if we are here, we have unused items in the kwargs dict, ie values for fields that don't exist
            raise ValueError("field names provided to 'add_row' for field names not contained in the update_fields list")
        
        # add a row of update data with the id_field value last in the list
        self.field_data.append( row_data )
        self.id_field_values.append( id_value )
        n_updated = 0
        if len( self.field_data ) >= self.buffer_size:
            n_updated = self._write_buffer()

        return n_updated
        
    def _write_buffer(self):
        """
        Commits the buffer to the database and empties it out. Called by add_row whenever we hit the 
        max buffer size of rows and by __exit__ to empty out the remaining buffer before we lose the 
        context.  
        Returns number of rows inserted
        """
        n = len(self.field_data)
        self._db.bulk_update(self.table_name, self.update_fields, self.id_field, self.field_data, self.id_field_values )
        
        # reset the buffer
        self.field_data = []
        self.id_field_values = []
        return n