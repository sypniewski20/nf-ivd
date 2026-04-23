#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ============================================================
// CLINICAL IVD GERMLINE PIPELINE
// ============================================================

include { Read_samplesheet; Read_bam }      from './modules/functions.nf'
include { dragmap_workflow }      from './subworkflows/mapping.nf' 
include { fastq_QC_workflow; mosdepth_workflow }     from './subworkflows/qc.nf'
include { multiqc_workflow }      from './subworkflows/multiqc.nf'
include { hc_workflow }           from './subworkflows/HC.nf'
include { manta_workflow }        from './subworkflows/manta.nf'
include { gCNV_workflow }         from './subworkflows/gCNV.nf'
include { germline_calibration_workflow }  from './subworkflows/calibration.nf'

def run_modes = params.run_mode?.split(',')*.trim()

workflow {

// 1. INPUT LAYER
    ch_fastqc_reports = Channel.empty()
    ch_mapping_stats  = Channel.empty()
    ch_mosdepth       = Channel.empty()
    ch_bam            = Channel.empty()

    if (params.input_type == 'fastq') {
        
        // --- NEW LOGIC FOR CALIBRATION / STREAMING MODE ---
        if ('calibration' in run_modes) {
            // 1. Parse the NIST-style samplesheet (URLs + MD5 as RGID)
            ch_raw_stream = Read_samplesheet(params.samplesheet)

            // 2. Stream through FASTP (Input: URLs -> Output: Local Filtered Chunks)

            // 4. Align grouped chunks into single Sample BAMs
            mapping_results = dragmap_workflow(ch_raw_stream)

        } else {
            // --- STANDARD LOCAL FASTQ MODE ---
            ch_fq = Read_samplesheet(params.samplesheet)
            ch_fq_qc = fastq_QC_workflow(ch_fq)

            filtered_fq = ch_fq_qc.fastq
            qc_results = ch_fq_qc.fastp
            ch_fastqc_reports = ch_fq_qc.fastqc

            mapping_results = dragmap_workflow(filtered_fq)
        }

        // --- COMMON POST-MAPPING LAYER ---
        mosdepth_results = mosdepth_workflow(mapping_results.ch_bam)
        
        ch_bam           = mapping_results.ch_bam
        ch_mapping_stats = mapping_results.ch_flagstat.collect()
        ch_mosdepth      = mosdepth_results.ch_mosdepth.collect()

    } else if (params.input_type == 'bam') {
        ch_bam = Read_bam(params.samplesheet)
    }

    // 2. ROUTING LAYER

    ch_final_vcf = Channel.empty()
    ch_final_stats = Channel.empty()

    hc_results = null

    if ('HC' in run_modes || 'calibration' in run_modes) {
        
        hc_results = hc_workflow(ch_bam)
        ch_final_vcf = hc_results.hc_vcf
        ch_vcf_stats = hc_results.stats
        
    }

    if ('SV' in run_modes) {
        manta_workflow(ch_bam)
        gCNV_workflow(ch_bam)

    }

    if ('calibration' in run_modes) {

        germline_calibration_workflow(hc_results.hc_vcf)

    }

    // 3. FINAL QC & REPORTING

    multiqc_workflow(
        ch_fastqc_reports.ifEmpty([]),
        ch_mapping_stats.ifEmpty([]), 
        ch_mosdepth.ifEmpty([]),
        ch_vcf_stats.ifEmpty([])
    )
}