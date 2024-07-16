# kraken-flu

- [kraken-flu](#kraken-flu)
  - [Purpose of this tool](#purpose-of-this-tool)
  - [Installation](#installation)
  - [Example Usage](#example-usage)
    - [Download the NCBI taxonomy files](#download-the-ncbi-taxonomy-files)
    - [Download sequence data](#download-sequence-data)
    - [Run the kraken-flu tool](#run-the-kraken-flu-tool)
    - [Use the kraken-flu outputs to create a kraken2 DB](#use-the-kraken-flu-outputs-to-create-a-kraken2-db)
  - [Detailed usage recipes](#detailed-usage-recipes)
  - [How it works](#how-it-works)
    - [Modifications to the influenza taxonomy](#modifications-to-the-influenza-taxonomy)
    - [Overview of the inner workings](#overview-of-the-inner-workings)
      - [Influenza genome completeness filter](#influenza-genome-completeness-filter)
    - [Note on RAM usage](#note-on-ram-usage)

## Purpose of this tool
This tool performs pre-processing of viral reference genome data for building custom kraken2 databases. It mostly deals with Influenza genomes, for which we create custom levels of hierarchy in the taxonomy.

The main tasks of the tool are:
- remove incomplete influenza genomes
- impose a unified naming scheme for influenza genomes
- re-organise the taxonomy for influenza genomes to create artificial taxa for segments (treating segments like separate species)

For details, see [How it works](#how-it-works).

## Installation
Install with pip:
```shell
pip install kraken_flu@git+ssh://git@gitlab.internal.sanger.ac.uk/malariagen1/misc_utils/kraken_flu.git
```

## Example Usage
### Download the NCBI taxonomy files
This can be done using the kraken2 tool kraken2-build like this (not part of this tool):
```shell
    kraken2-build --download-taxonomy --db TAXONOMY_PATH
```

### Download sequence data
As an example, download the viral section of the NCBI RefSeq resource into a directory SEQ_PATH:
```shell
wget -P [SEQ_PATH] https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
gunzip [SEQ_PATH]/viral.1.1.genomic.fna.gz
```

### Run the kraken-flu tool
This command creates the new taxonomy and sequence files in a temporary directory from which the kraken2-build tool can add them to the new database.  
The parameters are:  
--taxonomy_path: path to the NCBI taxonomy data on the local filesystem
--fasta_path: path to the FASTA file of genome sequences
--out_dir: temporary output directory for the modified taxonomy and sequence files
--filter: apply the filter for incomplete flu genomes
--filter_except: do not apply filter to genomes with this string matching the header

```shell
kraken_flu \
    --taxonomy_path  TAXONOMY_PATH \
    --fasta_path SEQ_PATH/viral.1.1.genomic.fna.gz \
    --out_dir TMP_DIR \
    --filter
    --filter_except "A/Goose/Guangdong/1/96(H5N1)" \
```

For further help with kraken-flu:
```shell
kraken_flu -h
```

### Use the kraken-flu outputs to create a kraken2 DB
Create a folder DB_DIR for the new database and run the kraken2-build tool to add the manipulated sequence data into the database folder like this. The command refers to the file TMP_DIR/library/library.fna, which is the default naming of the kraken-flu output for the manipulated FASTA file.

```shell
kraken2-build \
    --add-to-library TMP_DIR/library/library.fna \
    --db DB_DIR
```

Add the NCBI accession ID to taxon ID file from the taxonomy download. This is a large file, it is recommended to use a symlink to the original file here such as:

```shell
ln -s TAXONOMY_PATH/taxonomy/nucl_gb.accession2taxid DB_DIRtaxonomy
```

Now run the kraken2 database build process:

```shell
kraken2-build \
    --build \
    --db DB_DIR
```

The database is now ready to use with the kraken2 command to map unknown reads to the taxonomy and create reports.

## Detailed usage recipes
See the included jupyter notebooks for a real-world recipe for building a kraken2 database with kraken_flu in the [jupyter_notebooks](jupyter_notebooks) folder of this repo.

## How it works
The kraken_flu tool was created to support the GSU-ARD viral metagenomics pipeline by building custom kraken2 databases for simultaneous reference genome selection and typing/sub-typing of viruses. The focus of the tool is on Influenza taxonomy but it can also deal with other viruses that require modifications to the taxonomy in order to be useful for our metagenomics pipeline.

### Modifications to the influenza taxonomy
The original un-modified taxonomy for Influenza A viruses from [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=36420&lvl=3&lin=f&keep=1&srchmode=1&unlock) looks like this (showing two sub-types H1N1 and H2N3):
```
                          ┌─────────────┐                                   
                       ┌──┤Influenza A  ├──┐                                
                       │  └─────────────┘  │                                
                       │                   │                                
          ┌────────────▼─┐          ┌──────▼───────┐                        
          │ H1N1 subtype │          │ H3N2 subtype │ [...]                  
          └─┬────────────┘          └───────┬──────┘                        
            │                               │                               
            │ ┌──────────────────────────┐  │ ┌──────────────────────────┐  
            ├─► isolate genome sequence  │  ├─► isolate genome sequence  │  
            │ └──────────────────────────┘  │ └──────────────────────────┘  
            │                               │                               
            │ ┌──────────────────────────┐  │ ┌──────────────────────────┐  
            └─► isolate genome sequence  │  └─► isolate genome sequence  │  
              └──────────────────────────┘    └──────────────────────────┘  
                                                                            
                        [...]                            [...]              
```

In the ARD viral pipeline, we treat the 8 segments of the influenza virus like separate genomes, which enables the pipeline to detect reassortant viruses. Classifying segments instead of whole genomes also helps sub-typing influenza A because the sub-type defining segments 4 and 6 (HA and NA gene, respectively) are addressed in separation, taking away the noise from the other 6 segments, which do not define the subtype. 

The taxonomy imposed by kraken_flu looks like this for influenza A (not all nodes are shown):
```  
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
                                                                                 


```
Each genome is split into the 8 segments and each segment is assign as a child node of a new node that the tool adds to the flu taxonomy. 
For segments 4 and 6, the parental node is specific to the sub-type and segment (e.g. H1 segment 4, N2 segment 6 etc). All other segment sequences as well as the new sub-type/segment nodes for segments 4 and 6, are children of new taxonomy nodes for the virus type and segment number, e.g. "A segment 1".
The ARD viral pipeline uses the database created by this tool for the task of identifying the best-fitting reference sequence(s) for a sample and simultaneous typing/sub-typing. 

Take as an example a sample that contains an influenza A H1N1 virus. The viral pipeline uses the kraken2 tool with the database created by kraken_flu, which will show that the sample contains:
  -  influenza A segment 1, 2, 3, 5, 7 and 8
  -  influenza A H1 segment 4
  -  influenza A N1 segment 6

With this information, the pipeline can classify the virus in this sample as an influenza A H1N1 isolate. It would do so, even if, for example, segment 1 is the result of a reassortment event and stems from an H3N2 isolate, which could have skewed the result of a classification based on the whole genome.

### Overview of the inner workings
The two classes that do the heavy lifting are ```TaxonomyHandler``` and ```FastaHandler```. The process is orchestrated by a class ```KrakenDbBuilder``` and assisted by common functionality in ```utils```.
The ```TaxonomyHandler``` reads the NCBI taxonomy files ```nodes.dmp``` and ```names.dmp``` while ```FastaHanlder``` deals with the sequence files.  
The ```TaxonomyHandler``` class provides methods for adding the new taxonomy nodes as described in [Modifications to the influenza taxonomy](#modifications-to-the-influenza-taxonomy). It adds nodes for every possible segment of every possible type A and type B segment, regardless of whether we have any isolates to match them. This is in keeping with the overall design of kraken2 databases, which use the NCBI taxonomy - and therefore nodes for every organism classified in it - regardless of whether or not a particular database provides sequence data for each of those organisms.  
For influenza A, the new nodes include all possible sub-type nodes for segments 4 and 6. There are 18 known H subtypes and 11 known N subtypes. For influenza B, there are no sub-type specific segment parent nodes. Instead, all segment nodes for influenza B are at the same level.  
The ```TaxonomyHandler``` provides a method that identifies the correct new parent node for any influenza isolate segment by its name and attaches the isolate segment sequence to the correct new parent.
As an example: a sequence named "Influenza A (A/California/07/2009(H1N1)) segment 4" would be attached as a child node to the new taxonomy node "Influenza A H1 segment 4", in turn a child of "Influenza A segment 4".

__Caveat__: at present, the code relies on isolate names having an entry in the NCBI taxonomy names.dmp file. Thus, an isolate sequence in a FASTA file from outside of NCBI will not be added to the taxonomy right now. We are not currently using sequence sourced outside of NCBI.  
A solution to this problem is planned, see this [issue](https://gitlab.internal.sanger.ac.uk/malariagen1/misc_utils/kraken_flu/-/issues/8).

The ```FastaHandler``` class provides the methods for reading and writing FASTA files. The ```KrakenDbBuilder``` class uses a ```TaxonomyHandler``` and a ```FastaHandler``` object to create the modified taxonomy and cross-link the FASTA headers to the taxonomy before writing the modified taxonomy and sequence files to the output folder. 
Any sequences that are linked to new nodes in the taxonomy (currently only influenza) have a kraken:taxid tag added to the FASTA header before writing the output FASTA file. The resulting files serve as the input for the kraken2-build tool, which creates a kraken2 database. This tool uses the kraken:taxid tag, where present, to identify the matching node in the taxonomy tree. FASTA sequences that do not have a kraken:taxid tag are matched to the taxonomy tree using NCBI accession IDs. Thus, the output from kraken_flu for viruses that do not require modifications are not affected by passing through the tool and are added to the resulting kraken2 database as normal.

To run the tool, a CLI is provided by the ```cmd.py``` module. When installing the tool with pip, an executable command ```kraken_flu``` is installed locally. See [Run the kraken-flu tool](#run-the-kraken-flu-tool) for details on how to use it.

#### Influenza genome completeness filter
The final files are created by ```KrakenDbBuilder::create_db_ready_dir``` which has an option to run a filter for complete influenza genomes. 
If this filter is used, only those isolates will be kept, that have a full length copy of each of the 8 influenza genome segments. A list of isolate names can be provided, which are exempt from this filter. This is needed to keep, for example, the avian influenza reference Goose Guandong H5N1, which is incomplete.

### Note on RAM usage
For the sake of code simplicity and speed of execution, the input taxonomy and sequence files are read into RAM, which significantly simplifies the cross-referencing of the data between taxonomy and sequence library. The files involved are large and a run of this tool is memory hungry but this is not a limiting factor in the environment this was developed for. 
For reference: the taxonomy files from NCBI combined amount to around 400MB in RAM. The sequence libraries of viral RefSeq plus all Influenza genomes in Genbank amount to around 2GB. This should not pose a bottleneck in a modern data science environment.