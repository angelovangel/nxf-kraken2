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
params.database = "$HOME/db/minikraken_8GB_20200312"
params.help = ""

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
         --outdir           : ${params.outdir}
         --database         : ${params.database}

         Runtime data:
        -------------------------------------------
         Running with profile:   ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
         Running as user:        ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
         Launch dir:             ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
         Base dir:               ${ANSI_GREEN}${baseDir}${ANSI_RESET}
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
         --outdir           : where results will be saved, default is "results-fastp"
         --database         : kraken2 database, default is ${params.database}
        ===========================================
         """
         .stripIndent()

}
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
//println readcounts.size()

/* 
 * channel for kraken2 with reads, single- or pair-end 
 */

 Channel 
    .fromFilePairs( reads, checkIfExists: true, size: -1 ) // default is 2, so set to -1 to allow any number of files
    .ifEmpty { error "Can not find any reads matching ${reads}" }
    .set{ read_ch }

/* 
 * run kraken2 
 */

process kraken2 {

    //echo true
    publishDir params.outdir, mode: 'copy', pattern: 'fastp_trimmed/*' // publish 
    input:
        tuple sample_id, file(x) from read_ch
    
    output:
        file("kraken2.report")
        val seqmode into seqmode_ch


    script:
    def single = x instanceof Path // this is from Paolo: https://groups.google.com/forum/#!topic/nextflow/_ygESaTlCXg
    if ( !single ) {
        seqmode = "PE"
        """
        kraken2 -db ${params.database} \
        --fastq_input \
        --report kraken2.report
        --paired ${x[0]} ${x[1]} 
        """
    } 
    else {
        seqmode = "SE"
        """
        kraken2 -db ${params.database} \
        --fastq_input \
        --report kraken2.report
        ${x}
        """
    }

}
