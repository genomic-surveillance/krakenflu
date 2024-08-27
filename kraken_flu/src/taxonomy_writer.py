import os.path

from kraken_flu.src.db import Db

"""
This module contains the functionality for writing taxonomy data from the DB to files  
"""

def write_taxonomy(db: Db, output_dir:str):
    """
    Write the taxonomy files into output directory.  
    
    Args:
        output_dir: str, required
            Path to the outtput directory into which the nodes and names files are written.  
            The path is created if it doesn't exists yet.
    
    Returns:
        True on success
        
    Side-effects:
        Writes to file
            
    """
    if not os.path.exists( output_dir ):
        os.makedirs( output_dir )
        
    _write_nodes_file(db, os.path.join(output_dir,'nodes.dmp'))
    _write_names_file(db, os.path.join(output_dir,'names.dmp'))

    return True

def _write_nodes_file(db: Db, path:str):
    """
    Writes the nodes.dmp file.  
    Obtains an iterator over the database records from the Db object, so we are not 
    reading the large table into memory all at once.  

    Args:
        db: Db, required
            KrakenDbBuilder::Db object with the database connection
            
        path: str), required
            Path to write the nodes.dmp file
    """
    nodes_it = db.all_taxonomy_nodes_iterator()
    cols = [
        'tax_id', 
        'parent_tax_id', 
        'rank', 
        'embl_code', 
        'division_id', 
        'inherited_div_flag', 
        'genetic_code_id', 
        'inherited_GC_flag', 
        'mitochondrial_genetic_code_id', 
        'inherited_MGC_flag', 
        'GenBank_hidden_flag', 
        'hidden_subtree_root_flag', 
        'comments' 
    ]
    return _write_to_file_from_iterator(nodes_it, cols, path)
    
def _write_names_file(db: Db, path:str):
    """
    Writes the names.dmp file.  
    Obtains an iterator over the database records from the Db object, so we are not 
    reading the large table into memory all at once.  

    Args:
        db: Db, required
            KrakenDbBuilder::Db object with the database connection
            
        path: str), required
            Path to write the names.dmp file
    """
    names_it = db.all_taxonomy_names_iterator()
    cols = [
        'tax_id',
        'name',
        'unique_name',
        'name_class'
    ]
    return _write_to_file_from_iterator(names_it, cols, path)

def _write_to_file_from_iterator(it, cols, path):
    """
    Takes a database table results iterator and a file name and writes the data to the file.  
    The order of columns is dictated by the "cols" list of column/field names.  
    """
    with open( path, 'w' ) as out_fh:
        for item in it:
            rowdata = [ item[x] for x in cols]
            row_str = _format_for_tax_file_output( rowdata= rowdata)
            out_fh.write(row_str + "\n")
    return True

def _format_for_tax_file_output(rowdata:list):
    """
    Formats a list of data for a taxonomy file output row
    """
    return "\t|\t".join([ '' if x is None else str(x) for x in rowdata ]) + "\t|"