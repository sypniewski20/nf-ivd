include { DRAGMAP_BAM } from '../modules/mapping.nf'
include { MOSDEPTH; MOSDEPTH_EXOME } from '../modules/seqQC.nf'

// ==========================
// MAPPING WORKFLOW (PURE)
// ==========================
workflow mapping_workflow {

    take:
        ch_fastq
        fasta_dir
        fasta
        bed
        isWES   // <-- injected decision, no logic here

    main:

        DRAGMAP_BAM(ch_fastq, fasta_dir, fasta)

        bam      = DRAGMAP_BAM.out.ch_bam
        flagstat = DRAGMAP_BAM.out.ch_flagstat

        if (isWES) {
            MOSDEPTH_EXOME(bam, bed, fasta)
        } else {
            MOSDEPTH(bam, fasta)
        }

    emit:
        bam
        flagstat
        mosdepth = isWES ? MOSDEPTH_EXOME.out : MOSDEPTH.out
}