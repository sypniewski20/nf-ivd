#!/usr/bin/env nextflow

// ============================================================
// CLINICAL IVD GERMLINE PIPELINE
// Refactored for Nextflow >=26.04 (strict syntax parser, v2)
// ============================================================

include { readSamplesheet; readBam }                              from './modules/functions.nf'
include { dragmap_workflow }                                      from './subworkflows/mapping.nf'
include { fastq_QC_workflow; mosdepth_workflow; multiqc_workflow } from './subworkflows/qc.nf'
include { hc_workflow }                                            from './subworkflows/HC.nf'
include { deepvariant_workflow }                                   from './subworkflows/deepvariant.nf'
include { manta_workflow }                                         from './subworkflows/manta.nf'
include { delly_workflow }                                         from './subworkflows/delly.nf'
include { cnvkit_workflow }                                        from './subworkflows/cnvkit.nf'
include { gCNV_workflow }                                          from './subworkflows/gCNV.nf'
include { annotation_workflow }                                    from './subworkflows/annotation.nf'

workflow {

    // Was previously a bare top-level statement (outside any workflow/process
    // block). The strict parser, default since 26.04, does not allow script
    // declarations (include, workflow, ...) to be mixed with statements at
    // the top level -- statements must live inside a workflow. Also swapped
    // the spread-dot (*.trim()) for .collect{}, which is the pattern the
    // strict-syntax docs recommend, and named the closure parameter explicitly
    // instead of relying on implicit `it` (deprecated under the strict parser).
    def run_modes = params.run_mode?.split(',')?.collect { m -> m.trim() }

    // Clinical-safety fix: under the strict parser, CLI-supplied params are no
    // longer auto-cast to booleans/numbers. `--annotate false` now sets
    // params.annotate to the STRING 'false', which is truthy in Groovy, so
    // `params.annotate == true` or `if (params.annotate)` would silently run
    // annotation even when the user asked to disable it. Normalizing through
    // toString().toBoolean() handles both a real Boolean (from the config
    // file) and a CLI string override correctly.
    def annotate = params.annotate.toString().toBoolean()

    // 1. INPUT LAYER
    def ch_fastqc_reports = channel.empty()
    def ch_mapping_stats  = channel.empty()
    def ch_mosdepth       = channel.empty()
    def ch_bam            = channel.empty()

    if (params.input_type == 'fastq') {

        // --- STANDARD LOCAL FASTQ MODE ---
        def ch_fq    = readSamplesheet(params.samplesheet)
        def ch_fq_qc = fastq_QC_workflow(ch_fq)

        def filtered_fq = ch_fq_qc.fastq
        def qc_results   = ch_fq_qc.fastp
        ch_fastqc_reports = ch_fq_qc.fastqc

        def mapping_results = dragmap_workflow(filtered_fq)

        // --- COMMON POST-MAPPING LAYER ---
        def mosdepth_results = mosdepth_workflow(mapping_results.ch_bam)

        ch_bam           = mapping_results.ch_bam
        ch_mapping_stats = mapping_results.ch_metrics
        ch_mosdepth      = mosdepth_results.ch_mosdepth

        def multiqc_input = channel.empty()
            .mix(
                ch_mapping_stats,
                ch_mosdepth,
                qc_results
            )
            .flatten()
            .collect()

        multiqc_workflow(multiqc_input)

    } else if (params.input_type == 'bam') {
        ch_bam = readBam(params.samplesheet)

        ch_mosdepth = mosdepth_workflow(ch_bam).ch_mosdepth.collect()
        multiqc_workflow(ch_mosdepth)
    }

    def ch_final_vcf = channel.empty()
    def ch_final_tbi = channel.empty()

    def hc_results = null

    if ('HC' in run_modes) {

        hc_results = hc_workflow(ch_bam)

        if (annotate) {
            annotation_workflow(hc_results.ch_vcf, hc_results.ch_tbi)
        }

    }

    if ('DV' in run_modes) {

        def dv_results = deepvariant_workflow(ch_bam)

        if (annotate) {
            annotation_workflow(dv_results.ch_vcf, dv_results.ch_tbi)
        }

    }

    if ('SV' in run_modes) {
        manta_workflow(ch_bam)
        delly_workflow(ch_bam)
        cnvkit_workflow(ch_bam)
    }

    if ('gCNV' in run_modes) {
        gCNV_workflow(ch_bam)
    }

}
