include { DRAGMAP_BAM } from '../modules/mapping.nf'

// ==========================
// MULTIQC
// ==========================
workflow dragmap_workflow {

    take:
        ch_fq
    main:

        ch_fasta = Channel.value([
        file(params.fasta).parent,
        file(params.fasta),
        file("${params.fasta}.fai")
        ])

        DRAGMAP_BAM(ch_fq, ch_fasta)
    
    emit:

        // standardised outputs
        ch_bam       = DRAGMAP_BAM.out.ch_bam
        ch_flagstat  = DRAGMAP_BAM.out.ch_flagstat
        ch_md5       = DRAGMAP_BAM.out.ch_md5
}