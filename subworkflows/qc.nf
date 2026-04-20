include { FASTP_PROCESSING; FASTP_STREAM; FASTQC; MOSDEPTH; MOSDEPTH_EXOME } from '../modules/seqQC.nf'

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

workflow nist_streaming_QC_workflow {
    take:
        ch_nist_urls

    main:
        FASTP_STREAM(ch_nist_urls)
        
    emit:
        fastq  = FASTP_STREAM.out.fastq_filtered
        fastp  = FASTP_STREAM.out.fastp_log
}

workflow mosdepth_workflow {
    take:
        ch_bam

    main:
        def isWES = (params.seq_type == 'WES')

        ch_fasta = Channel.value([
        file(params.fasta),
        file("${params.fasta}.fai")
        ])

        if (isWES) {
            MOSDEPTH_EXOME(ch_bam, file(params.bed), ch_fasta) | set { ch_mosdepth }
        } else {
            MOSDEPTH(ch_bam, ch_fasta) | set { ch_mosdepth }
        }
    
    emit:

        // standardised outputs
        ch_mosdepth
}