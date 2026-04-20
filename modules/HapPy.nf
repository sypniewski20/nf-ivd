// ============================================================
// HAP.PY BENCHMARKING MODULES (CLINICAL IVD)
// ============================================================

process HAPY_GERMLINE_EVAL {
    label 'xxlarge'
    tag "${sample}"
    label "happy"

    publishDir "${params.outfolder}/${params.runID}/benchmark/happy/${sample}",
        mode: 'copy',
        overwrite: true

    input:
        // Input from combine: [sample_id, query_vcf, query_tbi, truth_vcf, truth_bed]
        tuple val(sample), path(query_vcf), path(query_tbi), path(truth_vcf), path(truth_bed)
        tuple path(fasta), path(fai)

    output:
        tuple val(sample), path("happy_${sample}"), emit: folder
        path("happy_${sample}/confident.summary.csv"), emit: summary_csv
        path("happy_${sample}/confident.vcf.gz"), emit: combined_vcf

    script:
    """
    set -euo pipefail
    mkdir -p happy_${sample}/global
    mkdir -p happy_${sample}/confident

    # 1. Global GIAB evaluation (Full Genome)
    happy \
        ${truth_vcf} \
        ${query_vcf} \
        -r ${fasta} \
        -o happy_${sample}/global/run \
        --engine vcfeval \
        --threads ${task.cpus} \
        --no-roc

    # 2. Confident region constrained evaluation (The 'Real' Metric)
    happy \
        ${truth_vcf} \
        ${query_vcf} \
        -f ${truth_bed} \
        -r ${fasta} \
        -o happy_${sample}/confident/run \
        --engine vcfeval \
        --threads ${task.cpus} \
        --no-roc
        
    # Standardize output names for the nextflow aggregator
    mv happy_${sample}/confident/run.summary.csv happy_${sample}/confident.summary.csv
    mv happy_${sample}/confident/run.vcf.gz happy_${sample}/confident.vcf.gz
    """
}

process HAPY_STRATIFIED_EVAL {
    label 'xxlarge'
    // Now we use the clean name passed from the workflow
    tag "${sample}:${strat_name}"
    label "happy"

    publishDir "${params.outfolder}/${params.runID}/benchmark/happy/${sample}/stratified/${strat_name}",
        mode: 'copy',
        overwrite: true

    input:
        // We added val(strat_name) here
        tuple val(sample), path(query_vcf), path(query_tbi), path(truth_vcf), path(truth_bed), path(strat_bed), val(strat_name)
        tuple path(fasta), path(fai)

    output:
        tuple val(sample), path("${sample}_${strat_name}_stratified"), emit: folder

    script:
    """
    set -euo pipefail
    mkdir -p ${sample}_${strat_name}_stratified

    happy \
        ${truth_vcf} \
        ${query_vcf} \
        -f ${strat_bed} \
        -r ${fasta} \
        -o ${sample}_${strat_name}_stratified/run \
        --engine vcfeval \
        --threads ${task.cpus} \
        --no-roc
    """
}