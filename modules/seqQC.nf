process FASTP_PROCESSING {

    tag "${sample}"
    label 'core'
    label 'medium'

    publishDir(
        { "${params.outfolder}/${params.runID}/fastp/${sample}" },
        mode: 'copy',
        overwrite: true,
        pattern: 'fastp.*'
    )

    input:
    tuple val(sample), val(ID), val(LB), val(PL), val(PU), path(read_1), path(read_2)

    output:
    tuple val(sample), val(ID), val(LB), val(PL), val(PU),
          path("${sample}.filtered.R1.fq.gz"),
          path("${sample}.filtered.R2.fq.gz"),
          emit: fastq_filtered

    tuple path("${sample}_fastp.html"),
          path("${sample}_fastp.json"),
          emit: fastp_log

    script:
    """
    fastp \
        -i ${read_1} \
        -I ${read_2} \
        -o ${sample}.filtered.R1.fq.gz \
        -O ${sample}.filtered.R2.fq.gz \
        -w ${task.cpus} \
        --html ${sample}_fastp.html \
        --json ${sample}_fastp.json \
        --detect_adapter_for_pe
    """
}


process FASTP_STREAM {

    tag "${sample}_${ID}"
    label 'core'
    label 'tiny'

    shell '/bin/bash'

    publishDir(
        { "${params.outfolder}/${params.runID}/fastp/${sample}" },
        mode: 'copy',
        overwrite: true,
        pattern: 'fastp.*'
    )

    input:
    tuple val(sample), val(ID), val(LB), val(PL), val(PU),
          path(R1_URL), path(R2_URL)

    output:
    tuple val(sample), val(ID), val(LB), val(PL), val(PU),
          path("${sample}_${ID}_R1_fastp.fq.gz"),
          path("${sample}_${ID}_R2_fastp.fq.gz"),
          emit: fastq_filtered

    tuple path("${sample}_${ID}_fastp.html"),
          path("${sample}_${ID}_fastp.json"),
          emit: fastp_log

    script:
    """
    fastp \
        -i <(curl -sL "${R1_URL}") \
        -I <(curl -sL "${R2_URL}") \
        -o ${sample}_${ID}_R1_fastp.fq.gz \
        -O ${sample}_${ID}_R2_fastp.fq.gz \
        -w ${task.cpus} \
        --html ${sample}_${ID}_fastp.html \
        --json ${sample}_${ID}_fastp.json \
        --detect_adapter_for_pe
    """
}


process MOSDEPTH {

    tag "${sample}"
    label 'qc'
    label 'medium'

    publishDir(
        "${params.outfolder}/${params.runID}/BAMQC",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            "${sample}/${filename}"
        }
    )

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(fasta), path(fai)

    output:
    path('*')

    script:
    """
    mosdepth \
        -f ${fasta} \
        -n \
        --fast-mode \
        --by 500 \
        --threshold 1,10,20,30 \
        -t ${task.cpus} \
        ${sample} \
        ${bam}
    """
}


process MOSDEPTH_EXOME {

    tag "${sample}"
    label 'qc'
    label 'medium'

    publishDir(
        "${params.outfolder}/${params.runID}/BAMQC",
        mode: 'copy',
        overwrite: true
    )

    input:
    tuple val(sample), path(bam), path(bai)
    path(intervals)
    tuple path(fasta), path(fai)

    output:
    path('*')

    script:
    """
    mosdepth \
        -f ${fasta} \
        -n \
        --fast-mode \
        --by ${intervals} \
        --threshold 1,10,20,30 \
        -t ${task.cpus} \
        ${sample} \
        ${bam}
    """
}


process FASTQC {

    tag "${sample}"
    label 'qc'
    label 'small'

    publishDir(
        { "${params.outfolder}/${params.runID}/fastqc/${sample}" },
        mode: 'copy',
        overwrite: true
    )

    input:
    tuple val(sample), val(ID), val(LB), val(PL), val(PU),
          path(read_1), path(read_2)

    output:
    path('*')

    script:
    """
    fastqc \
        -t ${task.cpus} \
        ${read_1} \
        ${read_2}
    """
}


process MULTIQC {

    label 'qc'
    label 'tiny'

    publishDir(
        "${params.outfolder}/${params.runID}/multiqc",
        mode: 'copy',
        overwrite: true
    )

    input:
    path(input)

    script:
    """
    multiqc . \
        --filename ${params.runID} \
        --verbose
    """
}