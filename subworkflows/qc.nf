include { FASTP_PROCESSING; FASTQC } from '../modules/seqQC.nf'

// ==========================
// QC WORKFLOW
// ==========================
workflow fastq_QC_workflow {

    take:
        ch_fq

    main:

        FASTP_PROCESSING(ch_fq)
        FASTQC(ch_fq)

    emit:

        // standardised outputs
        fastq   = FASTP_PROCESSING.out.fastq_filtered
        fastp   = FASTP_PROCESSING.out.fastp_log
        fastqc  = FASTQC.out
}