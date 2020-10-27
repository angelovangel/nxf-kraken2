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
params.fqpattern = "*_R{1,2}.fastq.gz"
params.readlen = 150
params.ontreads = false
params.kraken_db = "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20200919.tar.gz"
params.kraken_store = "$HOME/db/kraken"
params.kaiju_db = false
params.weakmem = false
params.taxlevel = "S" //level to estimate abundance at [options: D,P,C,O,F,G,S] (default: S)
params.skip_krona = false
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
         --readlen          : ${params.readlen}
         --outdir           : ${params.outdir}
         --kraken_db        : ${params.kraken_db}
         --kaiju_db         : ${params.kaiju_db}
         --weakmem          : ${params.weakmem}
         --taxlevel         : ${params.taxlevel}
         --skip_krona       : ${params.skip_krona}

         Runtime data:
        -------------------------------------------
         Running with profile:   ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
         Container:              ${ANSI_GREEN}${workflow.container}${ANSI_RESET}
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
         --readsdir     : directory with fastq files, default is "fastq"
         --fqpattern    : regex pattern to match fastq files, default is "*_R{1,2}.fastq.gz"
         --ontreads     : logical, set to true in case of Nanopore reads, default is false. This parameter has influence on fastp -q and bracken -r
         --readlen      : read length used for bracken, default is 150 (250 if ontreads is true). A kmer distribution file for this length has to be present in your database, see bracken help.
         --outdir       : where results will be saved, default is "results-kraken2"
         --kraken_db    : absolute path or ftp:// of kraken2 database, default is ${params.kraken_db}. See https://benlangmead.github.io/aws-indexes/k2 for available indexes
         --kaiju_db     : either 'false' (default, do not execute kaiju), or one of 'refseq', 'progenomes', 'viruses', 'nr' ...
         --weakmem      : logical, set to true to avoid loading the kraken2 database in RAM (on weak machines)
         --taxlevel     : taxonomical level to estimate bracken abundance at [options: D,P,C,O,F,G,S] (default: S)
        ===========================================
         """
         .stripIndent()

}

// because all is from conda, get versions from there. Note double escape for OR
process SoftwareVersions {
    publishDir "${params.outdir}/software_versions", mode: 'copy'

    output:
        file("software_versions.txt")

    script:
    """
    echo "software\tversion\tbuild\tchannel" > tempfile
    
    conda list | \
    grep 'fastp\\|kraken2\\|bracken\\|krona\\|r-data.table\\|r-dplyr\\|r-tidyr\\|r-dt\\|r-d3heatmap\\|r-base' \
    >> tempfile

    echo "kaiju \$(cd /kaiju && git tag | tail -n 1)" >> tempfile
    echo 'nextflow\t${nextflow.version}\t${nextflow.build}' >> tempfile
    multiqc --version | sed 's/, version//' >> tempfile

    # replace blanks with tab for easier processing downstream
    tr -s '[:blank:]' '\t' < tempfile > software_versions.txt
    """
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
process Fastp {

    tag "fastp on $sample_id"
    //echo true
    publishDir "${params.outdir}/trimmed_fastq", mode: 'copy', pattern: 'trim_*' // publish only trimmed fastq files

    input:
        tuple sample_id, file(x) from read_ch
    
    output:
        tuple sample_id, file('trim_*') into fastp_ch
        file("${sample_id}_fastp.json") into fastp4mqc_ch
        //tuple sample_id, file('trim_*') into fastp_ch


    script:
    def single = x instanceof Path // this is from Paolo: https://groups.google.com/forum/#!topic/nextflow/_ygESaTlCXg
    def fastp_input = single ? "-i \"${ x }\"" : "-i \"${ x[0] }\" -I \"${ x[1] }\""
    def fastp_output = single ? "-o \"trim_${ x }\"" : "-o \"trim_${ x[0] }\" -O \"trim_${ x[1] }\""
    def qscore_cutoff = params.ontreads ? 7 : 15 //here ontreads matters

    """
    fastp \
    -q $qscore_cutoff \
    $fastp_input \
    $fastp_output \
    -j ${sample_id}_fastp.json
    """
}
// make fastp channels for kraken2, kaiju and mqc
fastp_ch
    .into { fastp1; fastp2 }


if(params.kraken_db){
    Channel
        .of( "${params.kraken_db}" )
        .set { kraken_db }
} else {
        kraken_db = Channel.empty()
}

// setup kraken2 database, use url to tar.gz db
// input with ftp:// path downloads and stages the file
// input with absolute path stages the file downloaded previously
process KrakenDBPrep {
    storeDir "${params.kraken_store}"
    input:
        path kraken_file from kraken_db 
    output:
        path "**", type: 'dir' into kraken2_db_ch

// accomodate for cases where there is/there is not a leading directory in the tar archive
script:
"""
mkdir -p krakendb && tar -xf $kraken_file -C krakendb
"""
}


/* 
 run kraken2 AND bracken 
 Kraken-Style Bracken Report --> to use in pavian
 Bracken output file --> just a table, to be formatted and saved as html DataTable using R
 kraken2 counts file, this is the kraken2.output --> to use in krona
 */

process Kraken2 {
    tag "kraken2 on $sample_id"
    //echo true
    publishDir "${params.outdir}/samples", mode: 'copy', pattern: '*.{report,tsv}'
    
    input:
    // some channel work to accomodate for cases where there is/there is not a leading directory in the tar archive
    // see https://github.com/BenLangmead/aws-indexes/issues/5
        path db from kraken2_db_ch.flatten().last() 
        tuple sample_id, file(x) from fastp1
    
    output:
        file("*report") into kraken2mqc_ch // both kraken2 and the bracken-corrected reports are published and later used in pavian?
        file("*kraken2.krona") into kraken2krona_ch
        tuple sample_id, file("*bracken.tsv") into bracken2dt_ch
        file("*bracken.tsv") into bracken2summary_ch
    
    script:
    def single = x instanceof Path
    def kraken_input = single ? "\"${ x }\"" : "--paired \"${ x[0] }\"  \"${ x[1] }\""
    def memory = params.weakmem ? "--memory-mapping" : ""  // use --memory-mapping to avoid loading db in ram on weak systems
    def rlength = params.ontreads ? 250 : params.readlen // and here ontreads matters. Default for -r is 100 in bracken, Dilthey used 1k in his paper
    
        """
        kraken2 \
            -db $db \
            $memory \
            --report ${sample_id}_kraken2.report \
            $kraken_input \
            > kraken2.output
        cut -f 2,3 kraken2.output > ${sample_id}_kraken2.krona

        bracken \
            -d $db \
            -r $rlength \
            -i ${sample_id}_kraken2.report \
            -l ${params.taxlevel} \
            -o ${sample_id}_bracken.tsv
        """

}

// kaiju, in case params.kaiju_db is selected
// Conditionally create the input channel with data or as empty channel,
// the process consuming the input channels will only execute if the channel is populated
if(params.kaiju_db){
    Channel
        .of( "${params.kaiju_db}" )
        .set { kaiju_db }
} else {
        kaiju_db = Channel.empty()
}

process  KaijuDBPrep {
  input:
    val(x) from kaiju_db
  
  output:
    path("${x}/*.fmi") into fmi_ch
    path("*.dmp") into dmp_ch_1 // for kaiju
    path("*.dmp") into dmp_ch_2 // for kaiju2table
  
  script:
  """
  kaiju-makedb -s $x 
  """
}

// combine necessary here, for each sample_id with the db files 

// kaiju is also not executed if its input channels are empty
process Kaiju {
  //publishDir "${params.outdir}/samples", mode: 'copy', pattern: '*.tsv'

  input:
    tuple sample_id, file(x) from fastp2
    path("*") from dmp_ch_1.first() // this trick makes it a value channel, so no need to combine!
    file fmi from fmi_ch.first()
  
  output:
    file("*_kaiju.out") into kaiju_summary_ch
    file("*kaiju.out.krona") into kaiju2krona_ch

  script:
  def single = x instanceof Path
  def kaiju_input = single ? "-i \"${ x[0] }\"" : "-i \"${ x[0] }\" -j \"${ x[1] }\""
  """
  kaiju \
    -z 6 \
    -t nodes.dmp \
    -f $fmi \
    $kaiju_input \
    -o ${sample_id}_kaiju.out

  kaiju2krona -t nodes.dmp -n names.dmp -i ${sample_id}_kaiju.out -o ${sample_id}_kaiju.out.krona
  """
}

process KaijuSummary {
  publishDir params.outdir, mode: 'copy'
  
  input:
    path("*") from dmp_ch_2
    path x from kaiju_summary_ch.collect()
  
  output:
    path 'kaiju_summary.tsv' into kaiju2mqc_ch
  
  script:
  """
  kaiju2table \
    -t nodes.dmp \
    -n names.dmp \
    -r genus -m 1.0 \
    -o kaiju_summary.tsv \
    $x
  """
}

// setup the krona database and put it in a channel
// if no internet --skip_krona skips this process, and because the output channel is empty - it skips also
// process krona
    
process KronaDB {
    output:
        file("krona_db/taxonomy.tab") optional true into krona_db_ch // is this a value ch?

    when: 
        !params.skip_krona
        
    script:
    """
    ktUpdateTaxonomy.sh krona_db
    """
}

// prepare channel for krona, I want to have all samples in one krona plot
// e.g. ktImportTaxonomy file1 file2 ...

// run krona on the kraken2 and kaiju results
// SPLIT THIS PROCESS FOR KRAKEN AND KAIJU
process KronaFromKraken {
    publishDir params.outdir, mode: 'copy'

    input:
        file(x) from kraken2krona_ch.collect()
        //file(y) from kaiju2krona_ch.collect()
        file("krona_db/taxonomy.tab") from krona_db_ch
    
    output:
        file("*_taxonomy_krona.html")

    when:
        !params.skip_krona
    
    script:
    """
    mkdir krona
    ktImportTaxonomy -o kraken2_taxonomy_krona.html -tax krona_db $x
    """
}

process KronaFromKaiju {
    publishDir params.outdir, mode: 'copy'

    input:
        //file(x) from kraken2krona_ch.collect()
        file(y) from kaiju2krona_ch.collect()
        file("krona_db/taxonomy.tab") from krona_db_ch
    
    output:
        file("*_taxonomy_krona.html")

    when:
        !params.skip_krona
    
    script:
    """
    mkdir krona
    ktImportText -o kaiju_taxonomy_krona.html $y
    """
}

// format and save bracken table as DataTable, per sample

process DataTables1 {
    tag "DataTables1 on $sample_id"
    publishDir "${params.outdir}/samples", mode: 'copy', pattern: '*.html'

    input:
        tuple sample_id, file(x) from bracken2dt_ch
        
    output:
        file("*.html")

    script: 
    """
    bracken2dt.R $x ${sample_id}_bracken.html
    """
}

// use combine_bracken_outputs.py from bracken instead of bracken2summary.R
// should be the same though
// saves one summary table as html and csv for all samples
process DataTables2 {
    tag "DataTables2"
    publishDir params.outdir, mode: 'copy'

    input:
        file(x) from bracken2summary_ch.collect() //this gives all the bracken table files as params to the script
    output:
        file("*.html") optional true // for some sample numbers html is not generated
        file("*.csv")

// the bracken2summary.R decides what to output depending on number of samples
// for all sample counts - output a csv summary, for <= 12 - make html table and graph
    script:
    """
    bracken2summary.R $x
    """
}

// still no bracken module in mqc, opened an issue on github
process MultiQC {
    tag "MultiQC"
    publishDir params.outdir, mode: 'copy'

    input:
        file x from fastp4mqc_ch.collect()
        file y from kraken2mqc_ch.collect().ifEmpty([]) // do mqc if only one of kraken2 or kaiju was executed
        file z from kaiju2mqc_ch.ifEmpty([]) // do mqc if only one of kraken2 or kaiju was executed
    output:
        file "multiqc_report.html"
    
    script:
    """
    multiqc --interactive .
    """
}

//=============================
workflow.onComplete {
    if (workflow.success) {
        log.info """
            ===========================================
            Output files are here:   ==> ${ANSI_GREEN}$params.outdir${ANSI_RESET}
            ===========================================
            """
            .stripIndent()
    }
    else {
        log.info """
            ===========================================
            ${ANSI_RED}Finished with errors!${ANSI_RESET}
            """
            .stripIndent()
    }
}