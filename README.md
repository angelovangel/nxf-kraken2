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

Nothing to install, as soon as you have `docker` and `nextflow`. Download a `kraken2` database, e.g. from [here]() and run the pipeline:

```

```

## Download and setup of a `kraken2` database

