import os.path
import logging 
import subprocess

from kraken_flu.src.utils import parse_flu
from kraken_flu.src.db import Db

"""
This module contains the functionality for loading taxonomy data into the sqlite database    
"""

def load_taxonomy(db: Db, names_file_path:str, nodes_file_path:str, acc2taxid_file_path:str = None):
    """
    Main function. Orchestrates the loading of the NCBI taxonomy data into the database.

    Args:
        db: KrakenFlu::Db object, required
            This provides the connection to the DB
        
        names_file_path: str, required
            Path to the NCBI names.dmp file.
            
        nodes_path_file: str, required
            Path to the NCBI nodes.dmp file.  
            
        acc2taxid_file_path: str, optional
            Path to a NCBI accession 2 taxon ID file. This is optional but an assignment of all 
            sequences to taxon IDs in the FASTA header can only be done if this file is provided.
            
    Returns:
        True on success
        
    Side effects:
        Loads data into backend database
    """
    logging.info( f'starting to upload taxonomy data to DB')
    if not os.path.isfile( names_file_path ):
        raise ValueError(f"{names_file_path} is not a file")
    if not os.path.isfile( nodes_file_path ):
        raise ValueError(f"{nodes_file_path} is not a file")
    if acc2taxid_file_path and not os.path.isfile( acc2taxid_file_path ):
        raise ValueError(f"{acc2taxid_file_path} is not a file")
    
    _load_names(db,names_file_path)
    _load_nodes(db,nodes_file_path)
    if acc2taxid_file_path:
        _load_acc2taxids(db, acc2taxid_file_path)
    
    logging.info( f'Finished uploading taxonomy data to DB')
    return True

def _load_names(db:Db, names_file_path:str):
    """
    Upload the names.dmp file to the DB. Nodes and names are uploaded separately, relying on the NCBI 
    taxonomy file to link the two entities by tax_id, ie we are not checking that a taxonomy_node record exists for 
    the taxonomy_names we are inserting into the DB.
    Because this is a large number of records to insert, we are using a buffered insert method here.
    TODO: the only difference to _load_nodes is now the list of field names in the data dict.  Should combine the two methods.  
    """
    n_names = _get_num_records(names_file_path)
    logging.info( f'starting to upload {n_names} names records from {names_file_path} data to DB')
    
    with open( names_file_path, 'r' ) as fh:
        with db.bulk_insert_buffer(table_name='taxonomy_names', buffer_size= 50000) as b:
            for row in fh:
                d = _read_tax_data_file_row( row )
                n_inserted = b.add_row(
                    {
                        'tax_id': int(d[0]),
                        'name': d[1],
                        'unique_name': d[2],
                        'name_class': d[3]
                    }
                )
                if n_inserted > 0:
                    logging.info(f'flushed {n_inserted} records to DB')
            
    logging.info( f'finished uploading names records to DB')
    return True
    
def _load_nodes(db:Db, nodes_file_path:str):
    """
    Upload the names.dmp file to the DB
    """
    n_nodes = _get_num_records(nodes_file_path)
    logging.info( f'starting to upload {n_nodes} nodes records from {nodes_file_path} data to DB')
    
    with open( nodes_file_path, 'r' ) as fh:
        with db.bulk_insert_buffer(table_name='taxonomy_nodes', buffer_size= 50000) as b:
            for row in fh:
                d = _read_tax_data_file_row( row )
                n_inserted = b.add_row(
                    {
                        'tax_id': int(d[0]),
                        'parent_tax_id': int(d[1]),
                        'rank': d[2],
                        'embl_code': d[3],
                        'division_id': int(d[4]),                   
                        'inherited_div_flag': int(d[5]),            
                        'genetic_code_id': int(d[6]),
                        'inherited_GC_flag': int(d[7]),
                        'mitochondrial_genetic_code_id': int(d[8]),
                        'inherited_MGC_flag': int(d[9]),
                        'GenBank_hidden_flag': int(d[10]),
                        'hidden_subtree_root_flag': int(d[11]),
                        'comments': d[12]
                    }
                )
                if n_inserted > 0:
                    logging.info(f'flushed {n_inserted} records to DB')
            
    logging.info( f'finished uploading nodes records to DB')
    return True

def _load_acc2taxids(db:Db, acc2taxid_file_path:str):
    """
    Upload the accession to taxid file from NCBI to the DB
    """
    logging.info( f'starting to upload acc2taxid records from {acc2taxid_file_path} data to DB')
    raise NotImplemented("this function is not yet implemented and needs an update to the Db class: need a table for this data")
    logging.info( f'finish uploading {n} acc2taxid records to DB')
    
def _read_tax_data_file_row( row ):
    """
    Parses one row of data from names and nodes dmp file and returns as list
    Removes the trailing \t| from the last column
    """
    data = row.rstrip().split("\t|\t")
    data[-1] = data[-1].rstrip("\t|")
    return data

def _get_num_records( path ):
    """
    Get the number of records in the file by counting the lines.
    """
    p = subprocess.run(f"wc -l {path} | cut -d' ' -f1", shell=True, check=True, capture_output=True, encoding='utf-8')
    if not p.returncode:
        return int(p.stdout)
    else:
        raise Exception("failed to run wc -l to count lines in file")
