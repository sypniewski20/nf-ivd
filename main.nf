#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ============================================================
// CLINICAL IVD GERMLINE PIPELINE
// ============================================================

include { Read_samplesheet }      from './modules/functions.nf'
include { dragmap_workflow }      from './subworkflows/mapping.nf' 
include { fastq_QC_workflow; nist_streaming_QC_workflow; mosdepth_workflow }     from './subworkflows/qc.nf'
include { multiqc_workflow }      from './subworkflows/multiqc.nf'
include { hc_workflow }           from './subworkflows/HC.nf'
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
        ch_bam = Read_bam_checkpoint(params.samplesheet)
    }

    // 2. ROUTING LAYER 
    ch_final_vcf = Channel.empty()
    ch_final_stats = Channel.empty()

    switch(params.run_mode) {
        case 'HC':
            hc_results = hc_workflow(ch_bam)
            ch_final_vcf = hc_results.hc_vcf
            ch_final_stats = hc_results.stats

            break

        case 'calibration':
            
            hc_results = hc_workflow(ch_bam)

            cal_results = germline_calibration_workflow(hc_results.hc_vcf)
            // If calibration produces the benchmark VCF:
            ch_final_vcf = cal_results.vcf 
            ch_final_stats = cal_results.stats
            break

        default:
            error "CRITICAL: Unknown run_mode: ${params.run_mode}."
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