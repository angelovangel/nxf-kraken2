// kraken2-bracken-krona pipeline

/*
NXF ver 19.08+ needed because of the use of tuple instead of set
*/
if( !nextflow.version.matches('>=19.08') ) {
    println "This workflow requires Nextflow version 19.08 or greater and you are running version $nextflow.version"
    exit 1
}

/*
* ANSI escape codes to color output messages
*/
ANSI_GREEN = "\033[1;32m"
ANSI_RED = "\033[1;31m"
ANSI_RESET = "\033[0m"

/* 
 * pipeline input parameters 
 */
params.readsdir = "fastq"
params.outdir = "${params.readsdir}/results-kraken2" // output is where the reads are because it is easier to integrate with shiny later
params.fqpattern = "*_R{1,2}_001.fastq.gz"
params.ontreads = false
params.database = "$HOME/db/minikraken_8GB_20200312"
params.help = ""

/* 
 * handling of parameters 
 */

//just in case trailing slash in readsdir not provided...
readsdir_repaired = "${params.readsdir}".replaceFirst(/$/, "/") 
//println(readsdir_repaired)

// build search pattern for fastq files in input dir
reads = readsdir_repaired + params.fqpattern

// get counts of found fastq files
readcounts = file(reads)

if (params.help) {
    helpMessage()
    exit(0)
}

log.info """
        ===========================================
         K R A K E N 2 - B R A C K E N  P I P E L I N E

         Used parameters:
        -------------------------------------------
         --readsdir         : ${params.readsdir}
         --fqpattern        : ${params.fqpattern}
         --ontreads         : ${params.ontreads}
         --outdir           : ${params.outdir}
         --database         : ${params.database}

         Runtime data:
        -------------------------------------------
         Running with profile:   ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
         Running as user:        ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
         Launch dir:             ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
         Base dir:               ${ANSI_GREEN}${baseDir}${ANSI_RESET}
         Fastq files:            ${ANSI_GREEN}${ readcounts.size() } files found${ANSI_RESET}
         """
         .stripIndent()
/* 
 * define help 
 */
def helpMessage() {
log.info """
        ===========================================
         K R A K E N 2 - B R A C K E N  P I P E L I N E

         Note: 
         single- or pair-end data is automatically detected

         Usage:
        -------------------------------------------
         --readsdir         : directory with fastq files, default is "fastq"
         --fqpattern        : regex pattern to match fastq files, default is "*_R{1,2}_001.fastq.gz"
         --ontreads         : logical, set to true in case of Nanopore reads, default is false
         --outdir           : where results will be saved, default is "results-fastp"
         --database         : kraken2 database, default is ${params.database}
        ===========================================
         """
         .stripIndent()

}

/* 
 * channels for kraken2 with reads (single- or pair-end) and database
 */

Channel
    .fromFilePairs( reads, checkIfExists: true, size: -1 ) // default is 2, so set to -1 to allow any number of files
    .ifEmpty { error "Can not find any reads matching ${reads}" }
    .set{ read_ch }


/* 
 * run fastp 
 */
process fastp {

    tag "fastp on $sample_id"
    //echo true
    publishDir params.outdir, mode: 'copy', pattern: 'fastp_trimmed/*' // publish only trimmed fastq files

    input:
        tuple sample_id, file(x) from read_ch
    
    output:
        tuple sample_id, file('fastp_trimmed/trim_*') into fastp_ch


    script:
    def single = x instanceof Path // this is from Paolo: https://groups.google.com/forum/#!topic/nextflow/_ygESaTlCXg
    def qscore_cutoff = params.ontreads ? 7 : 15

    if ( !single ) {
        seqmode = "PE"
        """
        mkdir fastp_trimmed
        fastp \
        -q $qscore_cutoff \
        -i ${x[0]} -I ${x[1]} \
        -o fastp_trimmed/trim_${x[0]} -O fastp_trimmed/trim_${x[1]} \
        -j ${sample_id}_fastp.json
        """
    } 
    else {
        seqmode = "SE"
        """
        mkdir fastp_trimmed
        fastp \
        -q $qscore_cutoff \
        -i ${x} \
        -o fastp_trimmed/trim_${x} \
        -j ${sample_id}_fastp.json
        """
    }

}


/* 
 * run kraken2
 */
process kraken2 {
    tag "kraken2 on $sample_id"
    //echo true
    publishDir params.outdir, mode: 'copy', pattern: '*_kraken2.report'
    
    input:
        path krakendb from "${params.database}" //this db is not in the docker image
        tuple sample_id, file(x) from fastp_ch
    
    output:
        file("${sample_id}_kraken2.report") // this is published and later used in pavian?
        file("${sample_id}_kraken2.krona") into kraken2_ch
    
    script:
    def single = x instanceof Path

    if ( !single ) {
        """
        kraken2 \
            -db $krakendb \
            --report ${sample_id}_kraken2.report \
            --paired ${x[0]} ${x[1]} \
            > kraken2.output
        cut -f 2,3 kraken2.output > ${sample_id}_kraken2.krona
        """
    } 
    else {
        """
        kraken2 \
             -db ${params.database} \
            --report ${sample_id}_kraken2.report \
            ${x} \
            > kraken2.output
        cut -f 2,3 kraken2.output > ${sample_id}_kraken2.krona
        """
    }

}
// setup the krona database and put it in a channel
process krona_db {
    output:
        file("krona_db/taxonomy.tab") into krona_db_ch

    script:
    """
    ktUpdateTaxonomy.sh krona_db
    """
}

// prepare channel for krona, I want to have all samples in one krona plot
// e.g. ktImportTaxonomy file1 file2 ...

// run krona of the kraken2 results

process krona {
    publishDir params.outdir, mode: 'copy'

    input:
        file(x) from kraken2_ch.collect()
        file("krona_db/taxonomy.tab") from krona_db_ch
    
    output:
        file("*.html")
    
    script:
    """
    ktImportTaxonomy $x -o kraken2_taxonomy_krona.html -tax krona_db
    """
}