# nextflow-kraken2

A relatively simple metagenomics analysis pipeline written in nextflow. The pipeline is based on `kraken2` and `bracken` and is supplemented with `Krona` visualizations and interactive html tables.

## Description
The pipeline runs in a docker container by default. 
For a set of `fastq` files it executes:
- `fastp` 
- `kraken2` 
- `bracken` 
- `krona` plots are generated from the output of `bracken`

*Note: you have to setup your kraken2 database separately, it is not included in the container*

## Running the pipeline

Nothing to install, as soon as you have `docker` and `nextflow`. Download a `kraken2` database, e.g. from [here]() and run the pipeline:

```

```

## Download and setup of a `kraken2` database

