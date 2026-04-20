// ============================================================
// MANTA SUBWORKFLOWS
// ============================================================

include {
    JOINT_DIPLOID_ANALYSIS
} from '../modules/manta.nf'

workflow MANTA_WORKFLOW{
    take:
        ch_bam

    main:
        ch_fasta = Channel.value([
            file(params.fasta),
            file("${params.fasta}.fai")
            ])

        joint_bam = ch_bam.map{ it[0] }.collect()
        joint_bai = ch_bam.map{ it[1] }.collect()

        manta_results = JOINT_DIPLOID_ANALYSIS(joint_bam, joint_bai, ch_fasta)
    emit:
        manta_results
}