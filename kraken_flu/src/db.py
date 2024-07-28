from __future__ import annotations
from typing import Dict
from typing import Optional
from typing import List
from sqlalchemy import create_engine
from sqlalchemy import ForeignKey
from sqlalchemy import select
from sqlalchemy.orm import DeclarativeBase
from sqlalchemy.orm import Mapped
from sqlalchemy.orm import mapped_column
from sqlalchemy.orm import MappedAsDataclass
from sqlalchemy.orm import relationship
from sqlalchemy.orm import selectinload
from sqlalchemy.orm import Session
from sqlalchemy.orm.collections import attribute_keyed_dict
import os

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

class Base(DeclarativeBase):
    pass

# using the new Annotated Declarative Table style
# columns are declared with Python class on the left, mapped table column on the right
class TaxonomyNode(MappedAsDataclass, Base):
    __tablename__ = "taxonomy_nodes"

    # columns
    tax_id: Mapped[int] = mapped_column(primary_key=True, init=False)
    parent_tax_id: Mapped[Optional[int]] = mapped_column(
        ForeignKey("taxonomy_nodes.tax_id"), init=False
    )
    rank: Mapped[str]
    embl_code: Mapped[Optional[str]]
    division_id: Mapped[int]                   
    inherited_div_flag: Mapped[int]            
    genetic_code_id: Mapped[int]
    inherited_GC_flag: Mapped[int]
    mitochondrial_genetic_code_id: Mapped[int]
    inherited_MGC_flag: Mapped[int]
    GenBank_hidden_flag: Mapped[int]
    hidden_subtree_root_flag: Mapped[int]
    comments: Mapped[Optional[str]]

    # relationships
    taxonomy_names: Mapped[List[TaxonomyName]] = relationship(back_populates="taxonomy_node")
    sequence: Mapped[Optional[Sequence]] = relationship(back_populates="taxonomy_node")

    # children: delete-orphan will trigger a delete cascade when a parent is
    # deleted so that all children are deleted too
    children: Mapped[Dict[str, TaxonomyNode]] = relationship(
        cascade="all, delete-orphan",
        back_populates="parent",
        collection_class=attribute_keyed_dict("name"),
        init=False,
        repr=False,
    )

    parent: Mapped[Optional[TaxonomyNode]] = relationship(
        back_populates="children", remote_side=tax_id, default=None
    )
    
class TaxonomyName(MappedAsDataclass, Base):
    __tablename__ = "taxonomy_names"

# The code snippet you provided is defining the columns for the `TaxonomyNode` class using
# SQLAlchemy's declarative mapping. Let's break down the code snippet:
    tax_id: Mapped[Optional[int]] = mapped_column(
        ForeignKey("taxonomy_nodes.tax_id"), 
        primary_key=True,
        init=False
    )
    name: Mapped[str]= mapped_column(primary_key=True)
    name_class: Mapped[str]= mapped_column()
    unique_name: Mapped[Optional[str]] = mapped_column()

    # relationships
    taxonomy_node: Mapped[TaxonomyNode] = relationship(back_populates="taxonomy_names")


class Sequence(MappedAsDataclass, Base):
    __tablename__ = "sequences"

    # auto-increment ID primary key
    # init=False: A new object of this class does not require this attribute for its __init__
    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True, init=False)
    tax_id: Mapped[Optional[int]] = mapped_column(
        ForeignKey("taxonomy_nodes.tax_id"),
        init=False,
        nullable=True
    )
    fasta_header: Mapped[str]
    dna_sequence: Mapped[str]
    seq_length: Mapped[int]
    segment_number: Mapped[Optional[int]] = mapped_column(nullable=True)
    ncbi_acc: Mapped[Optional[str]] = mapped_column(nullable=True)
    flu_name: Mapped[Optional[str]] = mapped_column(nullable=True)
    flu_type: Mapped[Optional[str]] = mapped_column(nullable=True)
    flu_a_h_subtype: Mapped[Optional[int]] = mapped_column(nullable=True)
    flu_a_n_subtype: Mapped[Optional[int]] = mapped_column(nullable=True)
    include: Mapped[bool] = mapped_column()
    is_flu: Mapped[bool] = mapped_column()
    category: Mapped[str] = mapped_column(nullable=True)
    original_tax_id: Mapped[Optional[int]] = mapped_column(nullable=True)

    # relationships
    # 1:1 relationship to taxonomy_nodes: enforce single parent
    # NOTE: this uses a 1:1 relationship because the taxonomy and sequence data is loaded from different
    # files by different processes. Having all data in a single table would require linking the two types of
    # data at the moment we load them, which is not ideal. 
    taxonomy_node: Mapped[Optional[TaxonomyNode]] = relationship(back_populates="sequence", single_parent=True)


class Db():   
    """
    This class handles all interactions with the SQLite DB, including the DB session.

    Parameters:
        db_path: str, required
            absolute path to the sqlite file to be created
    """
    
    def __init__( self, db_path: str):
        self.db_path = db_path
        if os.path.exists(db_path):
            raise Exception(f"Cannot create database at `{db_path}` - the file exists already.")
        self._engine = create_engine(f"sqlite:///{db_path}", echo=False)
        # create the DB
        Base.metadata.create_all(self._engine)
        
        # create a session
        # This is a command line tool that does one job, so we use a single session throughout
        self._session = Session(self._engine)

    def add_sequence( self, fasta_header:str,  dna_sequence:str, category:str, flu_type:str, ncbi_acc:str, original_taxid:int, is_flu:bool, isolate_name:str, segment_number:int, h_subtype:int, n_subtype:int ):
        """
        Add a sequence record to table "sequences" without a link to a taxon node (which will be provided later).  
        NOTE on tax_id: If the FASTA header contains a kraken:taxid tag, we record this in the original_tax_id field.  
        The tax_id field is only populated later by this tool when we make the final decision about which node we are associating 
        this sequence with. So the original_tax_id is just a hint that may or may not end up being the tax_id we use in the end, but 
        both are recorded.
        """
        self._session.add(
            Sequence(
                fasta_header=fasta_header,
                dna_sequence=dna_sequence,
                seq_length=len(dna_sequence),
                ncbi_acc=ncbi_acc,
                is_flu=is_flu,
                flu_name=isolate_name,
                flu_type=flu_type,
                flu_a_h_subtype=h_subtype,
                flu_a_n_subtype=n_subtype,
                segment_number=segment_number,
                include=True,
                category=category,
                original_tax_id=original_taxid,
                taxonomy_node= None
            )
        )
        self._session.commit()
        
        
    def add_node(self, tax_id:int, parent_tax_id: int, rank: str, embl_code:str, division_id: int, inherited_div_flag:int, genetic_code_id:int,inherited_GC_flag:int, mitochondrial_genetic_code_id:int, inherited_MGC_flag:int, GenBank_hidden_flag:int, hidden_subtree_root_flag:int, comments:str):
        """
        Add a row into the nodes table.  
        """
        self._session.add(
            Node(
                tax_id= tax_id,
                parent_tax_id= parent_tax_id,
                rank= rank,
                embl_code= embl_code,
                division_id= division_id,                   
                inherited_div_flag= inherited_div_flag,            
                genetic_code_id= genetic_code_id,
                inherited_GC_flag= inherited_GC_flag,
                mitochondrial_genetic_code_id= mitochondrial_genetic_code_id,
                inherited_MGC_flag= inherited_MGC_flag,
                GenBank_hidden_flag= GenBank_hidden_flag,
                hidden_subtree_root_flag= hidden_subtree_root_flag,
                comments= comments
            )   
        )
        self._session.commit()

def add_name(self, ):
        """
        Add a row into the names table.  
        """
        self._session.add(
            Node(
                tax_id= tax_id,
                parent_tax_id= parent_tax_id,
                rank= rank,
                embl_code= embl_code,
                division_id= division_id,                   
                inherited_div_flag= inherited_div_flag,            
                genetic_code_id= genetic_code_id,
                inherited_GC_flag= inherited_GC_flag,
                mitochondrial_genetic_code_id= mitochondrial_genetic_code_id,
                inherited_MGC_flag= inherited_MGC_flag,
                GenBank_hidden_flag= GenBank_hidden_flag,
                hidden_subtree_root_flag= hidden_subtree_root_flag,
                comments= comments
            )   
        )
        self._session.commit()
