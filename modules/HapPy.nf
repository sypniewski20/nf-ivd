process HAPY_GERMLINE_EVAL {
    label 'xxlarge'
    tag "${sample}"
    label "happy"

    publishDir "${params.outfolder}/${params.runID}/benchmark/happy/${sample}",
        mode: 'copy',
        overwrite: true

    input:
        tuple val(sample), path(query_vcf)
        tuple path(truth_vcf), path(truth_bed)
        path(fasta)

    output:
        tuple val(sample), path("happy_${sample}")

    script:
    """
    set -euo pipefail

    mkdir -p happy_${sample}/global

    # ------------------------------------------------------------
    # Global GIAB evaluation
    # ------------------------------------------------------------
    happy \
        ${truth_vcf} \
        ${query_vcf} \
        -r ${fasta} \
        -o happy_${sample}/global \
        --engine vcfeval \
        --threads ${task.cpus} \
        --no-roc

    # ------------------------------------------------------------
    # Confident region constrained evaluation
    # ------------------------------------------------------------
    if [ -f ${truth_bed} ]; then
        happy \
            ${truth_vcf} \
            ${query_vcf} \
            -f ${truth_bed} \
            -r ${fasta} \
            -o happy_${sample}/confident \
            --engine vcfeval \
            --threads ${task.cpus} \
            --no-roc
    fi
    """
}

process HAPY_STRATIFIED_EVAL {
    label 'xxlarge'
    tag "${sample}:stratified"
    label "happy"

    publishDir "${params.outfolder}/${params.runID}/benchmark/happy/${sample}/stratified",
        mode: 'copy',
        overwrite: true

    input:
        tuple val(sample), path(query_vcf)
        tuple path(truth_vcf), path(truth_bed)
        path(strat_dir)
        path(fasta)

    output:
        path("${sample}_stratified")

    script:
    """
    set -euo pipefail

    mkdir -p ${sample}_stratified

    # ------------------------------------------------------------
    # Loop through all GIAB stratification BEDs
    # ------------------------------------------------------------
    for bed in ${strat_dir}/*.bed*; do

        strat=\$(basename \$bed .bed.gz)
        strat=\${strat%.bed}

        mkdir -p ${sample}_stratified/\${strat}

        happy \
            ${truth_vcf} \
            ${query_vcf} \
            -f \$bed \
            -r ${fasta} \
            -o ${sample}_stratified/\${strat} \
            --engine vcfeval \
            --threads ${task.cpus}

    done
    """
}