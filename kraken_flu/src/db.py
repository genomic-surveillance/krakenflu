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
    id: Mapped[int] = mapped_column(primary_key=True, init=False)
    tax_id: Mapped[Optional[int]] = mapped_column(
        ForeignKey("taxonomy_nodes.tax_id"), init=False
    )
    name: Mapped[str]= mapped_column()
    name_class: Mapped[str]= mapped_column()
    unique_name: Mapped[Optional[str]] = mapped_column()

    # relationships
    taxonomy_node: Mapped[TaxonomyNode] = relationship(back_populates="taxonomy_names")


class Sequence(MappedAsDataclass, Base):
    __tablename__ = "sequences"

    id: Mapped[int] = mapped_column(primary_key=True, init=False)
    tax_id: Mapped[Optional[int]] = mapped_column(
        ForeignKey("taxonomy_nodes.tax_id"), init=False
    )
    fasta_header: Mapped[str]
    dna_sequence: Mapped[str]
    seq_length: Mapped[int]
    ncbi_acc: Mapped[Optional[str]] = mapped_column()
    flu_name: Mapped[Optional[str]] = mapped_column()
    flu_a_h_subtype: Mapped[Optional[int]] = mapped_column()
    flu_a_n_subtype: Mapped[Optional[int]] = mapped_column()
    include: Mapped[bool] = mapped_column()

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
        os.makedirs( self.db_path )
        self._engine = create_engine("sqlite:///{db_path}", echo=False)
        # create the DB
        Base.metadata.create_all(self._engine)
        
        # create a session
        # This is a command line tool that does one job, so we use a single session throughout
        self._session = Session(self._engine)

        # more on adjacency list
        # https://docs.sqlalchemy.org/en/20/_modules/examples/adjacency_list/adjacency_list.html