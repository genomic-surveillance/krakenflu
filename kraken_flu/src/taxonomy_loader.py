import re
import os.path
import logging
import csv 

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
    if acc2taxid_file_path and not not os.path.isfile( acc2taxid_file_path ):
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
    """
    logging.info( f'starting to upload names from {names_file_path} data to DB')
    n=0
    with open( names_file_path, 'r' ) as fh:
        for row in fh:
            d = _read_tax_data_file_row( row )
            db.add_name(
                tax_id= int(d[0]),
                name= d[1],                
                unique_name= d[2],
                name_class = d[3]
            )
            n +=1            
    logging.info( f'finished uploading {n} names records to DB')
    return True
    
def _load_nodes(db:Db, nodes_file_path:str):
    """
    Upload the names.dmp file to the DB
    """
    logging.info( f'starting to upload names from {nodes_file_path} data to DB')
    n=0
    with open( nodes_file_path, 'r' ) as fh:
        for row in fh:
            d = _read_tax_data_file_row( row )
            db.add_node(
                tax_id= int(d[0]),
                parent_tax_id= int(d[1]),
                rank= d[2],
                embl_code= d[3],
                division_id= int(d[4]),                   
                inherited_div_flag= int(d[5]),            
                genetic_code_id= int(d[6]),
                inherited_GC_flag= int(d[7]),
                mitochondrial_genetic_code_id= int(d[8]),
                inherited_MGC_flag= int(d[9]),
                GenBank_hidden_flag= int(d[10]),
                hidden_subtree_root_flag= int(d[11]),
                comments= d[12]
            )
            n +=1
    logging.info( f'finished uploading {n} nodes records to DB')   
    return True

def _load_seq2taxids(db:Db, acc2taxid_file_path:str):
    """
    Upload the accession to taxid file from NCBI to the DB
    """
    logging.info( f'starting to upload names from {acc2taxid_file_path} data to DB')

    d_rename = {
        "accession":"ncbi_acc",
        "taxid":"tax_id"
    }

    n=0
    with open(acc2taxid_file_path, 'r') as fh:
        fd = csv.csv.DictReader(fh, delimiter="\t")(fh, delimiter="\t")
        for row in fd:
            row["version"] = row.pop("accession.version").split(".")[-1]
            for prior, post in d_rename.items():
                row[post] = row.pop(prior)

            db.add_seq2taxid(**row)
            
            n+=1
    logging.info( f'finish uploading {n} acc2taxid records to DB')
    return True
    
def _read_tax_data_file_row( row ):
    """
    Parses one row of data from names and nodes dmp file and returns as list
    Removes the trailing \t| from the last column
    """
    data = row.rstrip().split("\t|\t")
    data[-1] = data[-1].rstrip("\t|")
    return data