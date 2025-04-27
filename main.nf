#!/usr/bin/env nextflow

params.genomeFaPath = '/mainsd2/liujh/workspaces/data/common/hg19/hg19.fa'
params.genomeGtfPath = '/mainsd2/liujh/workspaces/data/common/hg19/genes/hg19.ncbiRefSeq.gtf'
params.sampleSheetPath = './test/params/sample_sheet_demo.tsv'
params.hisatIndexPrefix = '/mainsd2/liujh/workspaces/data/common/hg19/hisat/hg19'
params.hisat2IndexPrefix = '/mainsd2/liujh/workspaces/data/common/hg19/hisat2/hg19'
params.maxCpus = 10


workflow {
    rawCh = Channel
        .fromPath(params.sampleSheetPath)
        .splitCsv(header: true, sep: '\t', strip: true)
        .map { row ->
            tuple(
                [
                    id: row.id,
                    baseDirPath: rawGlobToBaseDirPath(row.rawFileGlob),
                    freecConfig: row.freecConfig,
                    species: row.species,
                ],
                file(row.rawFileGlob),
            )
        }

    QC_FASTP(rawCh)
    ALIGN_HISAT(QC_FASTP.out.reads, params.hisatIndexPrefix)
    QUANTIFY_FEATURE_COUNTS(ALIGN_HISAT.out.bam, params.genomeGtfPath)
}


process QC_FASTP {
    tag "${meta.id}"

    cpus { params.maxCpus > 16 ? 16 : params.maxCpus }

    memory '5 GB'

    storeDir "${meta.baseDirPath}/qc"
    fair true

    input:
    tuple val(meta), path(rawPath)

    output:
    tuple val(meta), path("${meta.id}_R{1,2}_qc.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}_fastp_report.html"), emit: html
    tuple val(meta), path("${meta.id}_fastp_report.json"), emit: json

    script:
    def String qcFastqPath1 = "${meta.id}_R1_qc.fastq.gz"
    def String qcFastqPath2 = "${meta.id}_R2_qc.fastq.gz"
    """
    #!/usr/bin/env bash
    fastp -i ${rawPath[0]} \
        -I ${rawPath[1]} \
        -o ${qcFastqPath1} \
        -O ${qcFastqPath2} \
        --html "${meta.id}_fastp_report.html" \
        --json "${meta.id}_fastp_report.json" \
        --thread ${task.cpus} \
        --detect_adapter_for_pe \
    """
}

process ALIGN_HISAT {
    tag "${meta.id}"
    cpus params.maxCpus
    memory '30 GB'
    // 30线程时，23GB
    storeDir "${meta.baseDirPath}/align"
    fair true

    input:
    tuple val(meta), path(qcFastqPath)
    val hisatIndexPrefix

    output:
    tuple val(meta), path("${meta.id}.bam", includeInputs: false), emit: bam
    tuple val(meta), path("${meta.id}_hisat_summary.txt", includeInputs: false), emit: summary

    script:
    def String bamPath = "${meta.id}.bam"
    def String hisatSummaryPath = "${meta.id}_hisat_summary.txt"
    """
    #!/usr/bin/env bash

    # hisat -S only support the output of SAM format.
    hisat -q \
        -x ${hisatIndexPrefix} \
        -1 ${qcFastqPath[0]} \
        -2 ${qcFastqPath[1]} \
        --threads ${task.cpus} \
        2>${hisatSummaryPath} |
        samtools view --threads ${task.cpus} -b -S - |
        samtools sort --threads ${task.cpus} -o ${bamPath} -
    # hisat -q \
    #     -x ${hisatIndexPrefix} \
    #     -1 ${qcFastqPath[0]} \
    #     -2 ${qcFastqPath[1]} \
    #     --threads ${task.cpus} \
    #     > ${bamPath} \
    #     2>&${hisatSummaryPath}
    """
}

process ALIGN_HISAT2 {
    tag "${meta.id}"
    cpus params.maxCpus
    storeDir "${meta.baseDirPath}/align_hisat2"
    fair true

    input:
    tuple val(meta), path(qcFastqPath)
    val hisat2IndexPrefix

    output:
    // TODO classify the outputs
    tuple val(meta), path("*", includeInputs: false)

    script:
    def String hisat2SummaryPath = "${meta.id}_hisat2_summary.bam"
    def String bamPath = "${meta.id}.bam"
    def String summaryPath = "${meta.id}_summary.txt"
    """
    #!/usr/bin/env bash

    hisat2 -q \
        -x ${hisat2IndexPrefix} \
        -1 ${qcFastqPath[0]} \
        -2 ${qcFastqPath[1]} \
        -S ${bamPath} \
        --summary-file ${summaryPath} \
        --threads ${task.cpus} \
        --mm
    """
}

process ALIGN_STAR {
    tag "${meta.id}"
    cpus params.maxCpus
    storeDir "${meta.baseDirPath}/align_hisat2"
    fair true

    input:
    tuple val(meta), path(qcFastqPath)
    val starGenomeDir

    output:
    // TODO classify the outputs
    tuple val(meta), path("*", includeInputs: false)

    script:
    """
    #!/usr/bin/env bash

    STAR --genomeDir ${starGenomeDir} \
        --readFilesIn ${qcFastqPath[0]} ${qcFastqPath[1]} \
        --outFileNamePrefix ./${meta.id} \
        --readFilesCommand gunzip -c \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped \
        --runThreadN ${task.cpus} \
        --genomeLoad LoadAndRemove
    """
}

process QUANTIFY_FEATURE_COUNTS {
    tag "${meta.id}"
    cpus { params.maxCpus > 8 ? 8 : params.maxCpus }
    memory '1 GB'
    storeDir "${meta.baseDirPath}/quantify"
    fair true

    input:
    tuple val(meta), path(bamPath)
    val genomeGtfPath

    output:
    tuple val(meta), path("${meta.id}.counts", includeInputs: false), emit: count
    tuple val(meta), path("${meta.id}.counts.summary", includeInputs: false), emit: summary

    script:
    def String countsPath = "${meta.id}.counts"
    """
    #!/usr/bin/env bash

    featureCounts -T ${task.cpus} \
        -F GTF \
        -p \
        -a ${genomeGtfPath} \
        -o ${countsPath} \
        ${bamPath}
    """
}

// 从原始文件匹配模式，计算得到对应的样本目录路径
def String rawGlobToBaseDirPath(rawFileGlob) {
    def String baseDirPath = rawFileGlob
    def String rawDirBasename = '/raw'
    if (!baseDirPath.contains(rawDirBasename)) {
        error("Invalid file pair glob: ${baseDirPath}")
    }
    return baseDirPath.substring(0, baseDirPath.lastIndexOf(rawDirBasename))
}
