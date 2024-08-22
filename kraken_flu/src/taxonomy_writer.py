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
            Path to the putput directory into which the nodes and names files are written
    
    Returns:
        True on success
        
    Side-effects:
        Writes to file
            
    """
    raise NotImplementedError

    return True

def _write_nodes_file(db: Db, path:str):
    """
    Writes the nodes.dmp file

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
    Writes the names.dmp file

    Args:
        db: Db, required
            KrakenDbBuilder::Db object with the database connection
            
        path: str), required
            Path to write the names.dmp file
    """
    names_it = db.all_taxonomy_names_iterator()
    cols = [
        'id', 
        'tax_id', 
        'name', 
        'name_class', 
        'unique_name'    
    ]
    return _write_to_file_from_iterator(names_it, cols, path)

def _write_to_file_from_iterator(it, cols, path):
    """
    Takes an iterator and a file name and writes the data to the file
    """
    with open( path, 'w' ) as out_fh:
        for node in it:
            rowdata = [ node[x] for x in cols]
            row_str = _format_for_tax_file_output( rowdata= rowdata)
            out_fh.write(row_str + "\n")
    return True

def _format_for_tax_file_output(rowdata:list):
    """
    Formats a list of data for a taxonomy file output row
    """
    return "\t|\t".join([ '' if x is None else str(x) for x in rowdata ]) + "\t|"