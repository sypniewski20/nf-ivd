include { DRAGMAP_BAM; DRAGMAP_STREAM_BAM; MERGE_BAM } from '../modules/mapping.nf'

def run_modes = params.run_mode?.split(',')*.trim()

workflow dragmap_workflow {

    take:
        ch_fq
    main:

        ch_fasta = Channel.value([
        file(params.fasta).parent,
        file(params.fasta),
        file("${params.fasta}.fai")
        ])

        if ('calibration' in run_modes) {

        ch_mapping = DRAGMAP_STREAM_BAM(ch_fq, ch_fasta)

        ch_grouped = ch_mapping
            .groupTuple()

        bams = ch_grouped.map { sample, lbs, bam, bai -> tuple(sample, lbs, bam) }
        bais = ch_grouped.map { sample, lbs, bam, bai -> tuple(sample, lbs, bai) }

        ch_bam_output = MERGE_BAM(bams, bais)

        } else {

            ch_bam_output = DRAGMAP_BAM(ch_fq, ch_fasta)

        }
    
    emit:

        // standardised outputs
        ch_bam       = ch_bam_output.ch_bam
        ch_flagstat  = ch_bam_output.ch_flagstat
        ch_md5       = ch_bam_output.ch_md5
}