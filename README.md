# nextflow-kraken2
a simple nextflow pipeline for running kraken2 and bracken in a docker container

## Description
The pipeline runs in a docker container by default. 
For a set of `fastq` files it executes:
- `fastp` 
- `kraken2` 
- `bracken` 
- `krona` plots are generated from the output of `bracken`

*Note: you have to setup your kraken2 database separately, it is not included in the container*

## Running the pipeline


## Setup `kraken2` database

