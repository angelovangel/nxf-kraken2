// kraken2-bracken-krona pipeline - DSL2

nextflow.enable.dsl = 2

/*
 * Check Nextflow version
 */
if (!nextflow.version.matches('>=21.04')) {
    println "This workflow requires Nextflow version 21.04 or greater -- You are running version $nextflow.version"
    exit 1
}

/*
 * ANSI escape codes to color output messages
 */
ANSI_GREEN = "\033[1;32m"
ANSI_RED = "\033[1;31m"
ANSI_RESET = "\033[0m"

/* 
 * Pipeline input parameters 
 */
params.readsdir = "fastq"
params.outdir = "${workflow.launchDir}/results-kraken2"
params.fqpattern = "*_R{1,2}.fastq.gz"
params.readlen = 150
params.ontreads = false
params.kraken_db = false
params.kaiju_db = false
params.weakmem = false
params.taxlevel = "S"
params.skip_krona = false
params.help = ""

/* 
 * Help message
 */
def helpMessage() {
    log.info """
        ===========================================
         K R A K E N 2 - B R A C K E N  P I P E L I N E

         Note: 
         Single- or paired-end data is automatically detected

         Usage:
        -------------------------------------------
         --readsdir     : directory with fastq files, default is "fastq"
         --fqpattern    : regex pattern to match fastq files, default is "*_R{1,2}.fastq.gz"
         --ontreads     : logical, set to true for Nanopore reads, default is false
         --readlen      : read length for bracken, default is 150 (250 if ontreads is true)
         --outdir       : where results will be saved, default is "results-kraken2"
         --kraken_db    : path to kraken2 database folder or 'false' to skip
         --kaiju_db     : 'refseq', 'progenomes', 'viruses', 'nr' or 'false' to skip
         --weakmem      : set to true to avoid loading kraken2 database in RAM
         --taxlevel     : taxonomical level for bracken [D,P,C,O,F,G,S] (default: S)
         --skip_krona   : skip making krona plots
        ===========================================
        """
        .stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// Repair readsdir path
readsdir_repaired = "${params.readsdir}".replaceFirst(/$/, "/")
reads = readsdir_repaired + params.fqpattern
readcounts = file(reads)

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
         Fastq files:            ${ANSI_GREEN}${readcounts.size()} files found${ANSI_RESET}
         """
         .stripIndent()

/*
 * PROCESSES
 */

process SoftwareVersions {
    publishDir "${params.outdir}/software_versions", mode: 'copy'

    output:
    path "software_versions.txt"

    script:
    """
    echo "software\tversion\tbuild\tchannel" > tempfile
    
    conda list | \
    grep 'fastp\\|kraken2\\|bracken\\|krona\\|r-data.table\\|r-dplyr\\|r-tidyr\\|r-dt\\|r-d3heatmap\\|r-base' \
    >> tempfile

    echo "kaiju \$(cd /kaiju && git tag | tail -n 1)" >> tempfile
    echo 'nextflow\t${nextflow.version}\t${nextflow.build}' >> tempfile
    multiqc --version | sed 's/, version//' >> tempfile

    tr -s '[:blank:]' '\t' < tempfile > software_versions.txt
    """
}

process Fastp {
    tag "fastp on $sample_id"
    publishDir "${params.outdir}/trimmed_fastq", mode: 'copy', pattern: 'trim_*'

    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path('trim_*'), emit: trimmed_reads
    path "${sample_id}_fastp.json", emit: json

    script:
    def single = reads instanceof Path
    def fastp_input = single ? "-i \"${reads}\"" : "-i \"${reads[0]}\" -I \"${reads[1]}\""
    def fastp_output = single ? "-o \"trim_${reads}\"" : "-o \"trim_${reads[0]}\" -O \"trim_${reads[1]}\""
    def qscore_cutoff = params.ontreads ? 7 : 15

    """
    fastp \
    -q $qscore_cutoff \
    $fastp_input \
    $fastp_output \
    -j ${sample_id}_fastp.json
    """
}

process Kraken2 {
    tag "kraken2 on $sample_id"
    publishDir "${params.outdir}/samples", mode: 'copy', pattern: '*.{report,tsv}'
    
    input:
    path db
    tuple val(sample_id), path(reads)
    
    output:
    path "*report", emit: reports
    path "*kraken2.krona", emit: krona
    tuple val(sample_id), path("*bracken.tsv"), emit: bracken_sample
    path "*bracken.tsv", emit: bracken_summary
    
    script:
    def single = reads instanceof Path
    def kraken_input = single ? "\"${reads}\"" : "--paired \"${reads[0]}\" \"${reads[1]}\""
    def memory = params.weakmem ? "--memory-mapping" : ""
    def rlength = params.ontreads ? 250 : params.readlen
    
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

process KaijuDBPrep {
    input:
    val db_name
  
    output:
    path "*.fmi", emit: fmi
    path "*.dmp", emit: dmp
  
    script:
    """
    kaiju-makedb -s $db_name
    """
}

process Kaiju {
    tag "kaiju on $sample_id"

    input:
    tuple val(sample_id), path(reads)
    path dmp
    path fmi
  
    output:
    path "*_kaiju.out", emit: kaiju_out
    path "*kaiju.out.krona", emit: krona

    script:
    def single = reads instanceof Path
    def kaiju_input = single ? "-i \"${reads[0]}\"" : "-i \"${reads[0]}\" -j \"${reads[1]}\""
    
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
    path dmp
    path kaiju_files
  
    output:
    path 'kaiju_summary.tsv', emit: summary
  
    script:
    """
    kaiju2table \
        -t nodes.dmp \
        -n names.dmp \
        -r genus -m 1.0 \
        -o kaiju_summary.tsv \
        $kaiju_files
    """
}

process KronaDB {
    output:
    path "krona_db", emit: db, optional: true

    when:
    !params.skip_krona
        
    script:
    """
    ktUpdateTaxonomy.sh krona_db
    """
}

process KronaFromKraken {
    publishDir params.outdir, mode: 'copy'

    input:
    path kraken_files
    path "krona_db"
    
    output:
    path "*_taxonomy_krona.html"

    when:
    !params.skip_krona
    
    script:
    """
    ktImportTaxonomy -o kraken2_taxonomy_krona.html -tax krona_db $kraken_files
    """
}

process KronaFromKaiju {
    publishDir params.outdir, mode: 'copy'

    input:
    path kaiju_files
    path "krona_db"
    
    output:
    path "*_taxonomy_krona.html"

    when:
    !params.skip_krona
    
    script:
    """
    ktImportText -o kaiju_taxonomy_krona.html $kaiju_files
    """
}

process DataTables1 {
    tag "DataTables1 on $sample_id"
    publishDir "${params.outdir}/samples", mode: 'copy', pattern: '*.html'

    input:
    tuple val(sample_id), path(bracken_file)
        
    output:
    path "*.html"

    script:
    """
    bracken2dt.R $bracken_file ${sample_id}_bracken.html
    """
}

process DataTables2 {
    tag "DataTables2"
    publishDir params.outdir, mode: 'copy'

    input:
    path bracken_files
    
    output:
    path "*.html", optional: true
    path "*.csv"

    script:
    """
    bracken2summary.R $bracken_files
    """
}

process MultiQC {
    tag "MultiQC"
    publishDir params.outdir, mode: 'copy'

    input:
    path fastp_json
    path kraken_reports
    path kaiju_summary
    
    output:
    path "multiqc_report.html"
    
    script:
    """
    multiqc --interactive .
    """
}

/*
 * WORKFLOW
 */

workflow {
    // Get software versions
    SoftwareVersions()
    
    // Create input channel from fastq files
    reads_ch = Channel
        .fromFilePairs(reads, checkIfExists: true, size: -1)
        .ifEmpty { error "Cannot find any reads matching ${reads}" }
    
    // Run Fastp
    Fastp(reads_ch)
    
    // Kraken2 workflow
    if (params.kraken_db) {
        kraken_db_ch = Channel.value(file(params.kraken_db))
        Kraken2(kraken_db_ch, Fastp.out.trimmed_reads)
        
        kraken_reports = Kraken2.out.reports
        kraken_krona = Kraken2.out.krona
        bracken_sample = Kraken2.out.bracken_sample
        bracken_summary = Kraken2.out.bracken_summary
        
        // DataTables for individual samples
        DataTables1(bracken_sample)
        
        // DataTables summary
        DataTables2(bracken_summary.collect())
    } else {
        kraken_reports = Channel.empty()
        kraken_krona = Channel.empty()
    }
    
    // Kaiju workflow
    if (params.kaiju_db) {
        kaiju_db_ch = Channel.value(params.kaiju_db)
        KaijuDBPrep(kaiju_db_ch)
        
        Kaiju(
            Fastp.out.trimmed_reads,
            KaijuDBPrep.out.dmp.collect(),
            KaijuDBPrep.out.fmi.first()
        )
        
        KaijuSummary(
            KaijuDBPrep.out.dmp.collect(),
            Kaiju.out.kaiju_out.collect()
        )
        
        kaiju_summary = KaijuSummary.out.summary
        kaiju_krona = Kaiju.out.krona
    } else {
        kaiju_summary = Channel.empty()
        kaiju_krona = Channel.empty()
    }
    
    // Krona plots
    if (!params.skip_krona) {
        KronaDB()
        
        if (params.kraken_db) {
            KronaFromKraken(
                kraken_krona.collect(),
                KronaDB.out.db
            )
        }
        
        if (params.kaiju_db) {
            KronaFromKaiju(
                kaiju_krona.collect(),
                KronaDB.out.db
            )
        }
    }
    
    // MultiQC
    MultiQC(
        Fastp.out.json.collect(),
        kraken_reports.collect().ifEmpty([]),
        kaiju_summary.ifEmpty([])
    )
}

/*
 * Completion handler
 */
workflow.onComplete {
    if (workflow.success) {
        log.info """
            ===========================================
            Output files are here:   ==> ${ANSI_GREEN}${params.outdir}${ANSI_RESET}
            ===========================================
            """
            .stripIndent()
    } else {
        log.info """
            ===========================================
            ${ANSI_RED}Finished with errors!${ANSI_RESET}
            """
            .stripIndent()
    }
}
