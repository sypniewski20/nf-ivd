
include { Read_samplesheet }       from './modules/functions.nf'
include { Read_bam_checkpoint }    from './modules/functions.nf'

include { preprocessing_workflow } from './subworkflows/preprocessing.nf'
include { fastq_QC_workflow }      from './subworkflows/qc.nf'
include { mapping_workflow }       from './subworkflows/mapping.nf'
include { multiqc }                from './subworkflows/multiqc.nf'

include { hc_workflow }     from './subworkflows/HC.nf'
include { tumor_only_workflow }   from './subworkflows/tumor_only.nf'
include { tumor_normal_workflow } from './subworkflows/tumor_normal.nf'
include { calibration_workflow }  from './subworkflows/calibration.nf'

// =========================
// INPUT LAYER
// =========================
if (params.input_type == 'fastq') {

    ch_fq = Read_samplesheet(params.samplesheet)

    qc = fastq_QC_workflow(ch_fq)

    ch_bam = mapping_workflow(
        qc.fastq,
        params.fasta_dir,
        params.fasta,
        params.bed
    )

} else if (params.input_type == 'bam') {

    ch_bam = Read_bam_checkpoint(params.samplesheet)
    qc = null

} else {
    error "Unknown input_type: ${params.input_type}"
}

// =========================
// ROUTING LAYER (RUN MODE)
// =========================
switch(params.run_mode) {

    case 'HC':
        result = hc_workflow(ch_bam, params)
        break

    case 'tumor_only':
        result = tumor_only_workflow(ch_bam, params)
        break

    case 'tumor_normal':
        result = tumor_normal_workflow(ch_bam, params)
        break

    case 'calibration':
        result = calibration_workflow(ch_bam, params)
        break

    default:
        error "Unknown run_mode: ${params.run_mode}"
}

// =========================
// FINAL QC
// =========================
if (qc != null) {
    multiqc(
        qc.fastqc,
        result.flagstat,
        result.mosdepth
    )
}