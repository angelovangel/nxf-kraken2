# nextflow-kraken2

A relatively simple metagenomics analysis pipeline written in nextflow [[1]](#1). The pipeline is based on `kraken2` and `bracken`, and is supplemented with `Krona` visualizations and interactive html tables.

## Description

The pipeline runs in a docker container by default. Both Illumina and Nanopore (but not hybrid) data can be processed. For a set of `fastq` files it executes:

- [`fastp`](https://github.com/OpenGene/fastp) - filter and trim reads with default parameters
- [`kraken2`](http://ccb.jhu.edu/software/kraken2/) [[2]](#2) - taxonomic assignment of the reads  (you have to setup your kraken2 database separately, it is not included in the container image)
- [`bracken`](http://ccb.jhu.edu/software/bracken/) [[3]](#3) - abundance estimation at a single level in the taxonomic tree, e.g. species
- [`krona`](https://github.com/marbl/Krona/wiki) [[4]](#4) - plots are generated from the output of `kraken2`
- [`DataTables`](https://datatables.net/) - generates an interactive HTML table with the results from `bracken` for each sample, as well as a summary table for all the samples

## Installation and running the pipeline

Nothing to install, as soon as you have `docker` and `nextflow`. Setup a `kraken2` database (see below), and run the pipeline:

```bash
# run with a test dataset (included)
nextflow run angelovangel/nextflow-kraken2 -profile test

# see options and how to run
nextflow run angelovangel/nextflow-kraken2 --help

```

For all the arguments and how to use them see the output of `nextflow run angelovangel/nextflow-kraken2 --help`

## Output

All output files are in the folder `results-kraken2`, which is in the folder with reads data used for running the pipeline. An example of the outputs, generated with a small Illumina dataset can be found under `example_output` in this repository.

The outputs are:

- `timmed_fastq/` - directory with fastq files after trimming, these are also used for taxonomic profiling
- `bracken_summary_heatmap/table.html`- standalone html files with summary information from bracken. Note that these files will be generated only if there are less than 24 samples
- `bracken_summary_long/wide.csv`- summary bracken information (all found taxa in all samples), in different formats
- `kraken2taxonomy_krona.html`-  an interactive Krona plot of the kraken2 output
- `samples/` - directory with individual (per sample) kraken2 and bracken-corrected report files (*.report) and with the abundance table from bracken (as html and tsv). *Tip: the report files can be directly imported in [Pavian](https://github.com/fbreitwieser/pavian) for nice interactive visualizations.

## Download and setup of a `kraken2` database

Get the minikraken2 database, e.g. from [here](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads), put it in a suitable folder ($HOME/db/) and untar:

```bash
tar -xzvf minikraken_8GB_202003.tgz
```

This pre-built minikraken2 database has the required Bracken files included (for read lengths 50, 100, 150, 200 and 250).

## References

This pipeline just uses some really nice work from others:   

<a id="1">[1]</a> 
P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017) https://doi.org/10.1038/nbt.3820

<a id="2">[2]</a> 
Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019) https://doi.org/10.1186/s13059-019-1891-0

<a id="3">[3]</a> 
Lu J, Breitwieser FP, Thielen P, Salzberg SL. 2017. Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science 3:e104 https://doi.org/10.7717/peerj-cs.104

<a id="4">[4]</a> 
Ondov BD, Bergman NH, Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011;12:385. Published 2011 Sep 30. doi:10.1186/1471-2105-12-385