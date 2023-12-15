# kraken_flu
This tool takes care of the modifications that need to be made to KRAKEN2 database input files in order to build a KRAKEN2 database where the segments of influenza genomes appear as new taxa, one per segment. 

The default way of handling such genomes in KRAKEN2 is one taxon per genome.

## Installation
Install with pip. You will probably want to create a venv for this first.
```shell
pip install kraken_flu@git+ssh://git@gitlab.internal.sanger.ac.uk/malariagen1/misc_utils/kraken_flu.git
```

## Usage
Once installed, run with the following command:

```shell
kraken_flu --taxonomy_path FILE --library_path FILE [--acc2tax_path ACC2TAX_PATH] --out_dir OUT_DIR
```

### Example database creation
Use kraken2 build tools to download taxonomy and genome sequence data to build a database. 

Step 1: obtain the taxonomy files. These are retrieved from [NCBI](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) by the kraken tool.
```shell
kraken2-build \
    --download-taxonomy \
    --db SOME/PATH/
```
This command creates a directory SOME/PATH/taxonomy

Step 2: download genome data. In this case, using kraken2 to download a prebuilt version of NCBI RefSeq for the virus division.  

```shell
kraken2-build \
    --download-library viral \
    --db SOME/PATH/
```
This command creates a directory SOME/PATH/library/viral

Step 3: Use the kraken_flu tool to scan the files created above and re-assign influenza H1N1 and H3N2 virus genome segments to new taxa, one segment per taxon ID.   
The tool will create a new directory with the modified versions of the taxonomy and library (FASTA) files. The paths are taken from the above example kraken2 build commands.  

```shell
kraken_flu \
    --taxonomy_path SOME/PATH/taxonomy \
    --library_path SOME/PATH/library/viral \
    --out_dir SOME/OUTPUT/DIR
```

Step 4: use kraken2 to build the database from the modified files  

```shell
kraken2-build \
    --build \
    --db  SOME/OUTPUT/DIR
```

Optional clean-up step: remove (potentially large) files that are no longer needed in the build directory.  

```shell
kraken2-build \
    --clean \
    --db SOME/OUTPUT/DIR
```