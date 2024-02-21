# kraken-flu
This tool performs pre-processing of viral reference genome data for building custom kraken2 databases.

The main tasks are:
- remove incomplete influenza genomes
- impose a unified naming scheme for influenza genomes
- re-organise the taxonomy for influenza genomes to create artificial taxa for segments (treating segments like separate species)

## Installation
Install with pip. You will probably want to create a venv for this first.
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
wget -P SEQ_PATH https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
gunzip SEQ_PATH/viral.1.1.genomic.fna.gz
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