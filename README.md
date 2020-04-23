# nextflow-kraken2

A relatively simple metagenomics analysis pipeline written in nextflow. The pipeline is based on `kraken2` and `bracken`, and is supplemented with `Krona` visualizations and interactive html tables.

## Description

The pipeline runs in a docker container by default. Both Illumina and Nanopore (but not hybrid) data can be processed. For a set of `fastq` files it executes:

- [`fastp`](https://github.com/OpenGene/fastp) - filter and trim reads with default parameters
- [`kraken2`](http://ccb.jhu.edu/software/kraken2/) - taxonomic assignment of the reads (you have to setup your kraken2 database separately, it is not included in the container image)
- [`bracken`](http://ccb.jhu.edu/software/bracken/) - abundance estimation at a single level in the taxonomic tree, e.g. species
- [`krona`](https://github.com/marbl/Krona/wiki) - plots are generated from the output of `kraken2`
- [`DataTables`](https://datatables.net/) - generates an interactive HTML table with the results from `bracken` for each sample, as well as a summary table for all the samples

## Running the pipeline

### Installation

Nothing to install, as soon as you have `docker` and `nextflow`. Setup a `kraken2` database (see below), and run the pipeline:

```bash
# run with a test dataset (included)
nextflow run angelovangel/nextflow-kraken2 -profile test

# see options and how to run
nextflow run angelovangel/nextflow-kraken2 --help

```

### Main arguments

See output of `nextflow run angelovangel/nextflow-kraken2 --help`

### Output

All output files are in the folder `results-kraken2`, which is in the folder with reads data used for running the pipeline. An example of the outputs, generated with Nanopore reads from the Loman's lab [Nanopore GridION Mock Microbial Community Data Community Release](https://github.com/LomanLab/mockcommunity) can be found under `example_output` in this repository.

## Download and setup of a `kraken2` database

Get the minikraken2 database, e.g. from [here](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads), put it in a suitable folder ($HOME/db/) and untar:

```bash
tar -xzvf minikraken_8GB_202003.tgz
```

The pre-built minikraken2 database has the required Bracken files included (for read lengths 50, 100, 150, 200 and 250).
