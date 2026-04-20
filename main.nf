#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ============================================================
// CLINICAL IVD GERMLINE PIPELINE
// ============================================================

include { Read_samplesheet; Read_bam }      from './modules/functions.nf'
include { dragmap_workflow }      from './subworkflows/mapping.nf' 
include { fastq_QC_workflow; nist_streaming_QC_workflow; mosdepth_workflow }     from './subworkflows/qc.nf'
include { multiqc_workflow }      from './subworkflows/multiqc.nf'
include { hc_workflow }           from './subworkflows/HC.nf'
include { manta_workflow }        from './subworkflows/manta.nf'
include { gCNV_workflow }         from './subworkflows/gCNV.nf'
include { germline_calibration_workflow }  from './subworkflows/calibration.nf'

workflow {

// 1. INPUT LAYER
    ch_fastqc_reports = Channel.empty()
    ch_mapping_stats  = Channel.empty()
    ch_mosdepth       = Channel.empty()
    ch_bam            = Channel.empty()

    if (params.input_type == 'fastq') {
        
        // --- NEW LOGIC FOR CALIBRATION / STREAMING MODE ---
        if (params.run_mode == 'calibration') {
            // 1. Parse the NIST-style samplesheet (URLs + MD5 as RGID)
            ch_raw_stream = Read_samplesheet(params.samplesheet)

            // 2. Stream through FASTP (Input: URLs -> Output: Local Filtered Chunks)
            nist_streaming_QC_workflow(ch_raw_stream)

            // 3. Group all MD5-tagged chunks by NIST_SAMPLE_NAME (HG002, etc.)
            // Output of nist_streaming_QC_workflow: [RGSM, RGID, RGLB, RGPL, R1_file, R2_file]
            ch_grouped_to_map = nist_streaming_QC_workflow.out.fastq
                .groupTuple(by: 0) 

            // 4. Align grouped chunks into single Sample BAMs
            mapping_results = dragmap_workflow(ch_grouped_to_map)
            
            // Collect QC from streaming fastp
            qc_results = nist_streaming_QC_workflow.out.fastp

        } else {
            // --- STANDARD LOCAL FASTQ MODE ---
            ch_fq = Read_samplesheet(params.samplesheet)
            ch_fq_qc = fastq_QC_workflow(ch_fq)

            filtered_fq = ch_fq_qc.fastq
            qc_results = ch_fq_qc.fastp
            fastqc_reports = ch_fq_qc.fastqc

            mapping_results = dragmap_workflow(filtered_fq)
        }

        // --- COMMON POST-MAPPING LAYER ---
        mosdepth_results = mosdepth_workflow(mapping_results.ch_bam)
        
        ch_bam           = mapping_results.ch_bam
        ch_mapping_stats = mapping_results.ch_flagstat
        ch_mosdepth      = mosdepth_results.ch_mosdepth

    } else if (params.input_type == 'bam') {
        ch_bam = Read_bam(params.samplesheet)
    }

    // 2. ROUTING LAYER

    ch_final_vcf = Channel.empty()
    ch_final_stats = Channel.empty()

    def run_modes = params.run_mode?.split(',')*.trim()

    hc_results = null

    if ('HC' in run_modes || 'calibration' in run_modes) {
        hc_results = hc_workflow(ch_bam)

        if ('HC' in run_modes) {
            ch_final_vcf = ch_final_vcf.mix(hc_results.hc_vcf)
            ch_final_stats = ch_final_stats.mix(hc_results.stats)
        }
    }

    if ('SV' in run_modes) {
        manta_results = manta_workflow(ch_bam)
        gcnv_results = gCNV_workflow(ch_bam)

        ch_final_vcf = ch_final_vcf.mix(manta_results.vcf)
        // optionally:
        // ch_final_vcf = ch_final_vcf.mix(gcnv_results.vcf)
    }

    if ('calibration' in run_modes) {

        cal_results = germline_calibration_workflow(hc_results.hc_vcf)

        ch_final_vcf = ch_final_vcf.mix(cal_results.vcf)
        ch_final_stats = ch_final_stats.mix(cal_results.stats)
    }

    // 3. FINAL QC & REPORTING
    // We mix the mapping stats with any tool-specific stats
    multiqc_workflow(
        ch_fastqc_reports.collect().ifEmpty([]),
        ch_mapping_stats.collect().ifEmpty([]), 
        ch_mosdepth.collect().ifEmpty([]),
        ch_final_stats.collect().ifEmpty([])
    )
}