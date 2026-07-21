include {
    VEP_GERMLINE_SNV; SPLICE_AI; FILTER_VEP
} from '../modules/annotations.nf'

workflow annotation_workflow {
    take:
        ch_vcf
        ch_tbi
    main:
        ch_fasta = Channel.value([
            file(params.fasta),
            file("${params.fasta}.fai")
            ])


        if (params.seq_type != "WES") {
            SPLICE_AI(ch_vcf, ch_tbi, ch_fasta)
            ch_spliceai_out = SPLICE_AI.out
        } else {
            ch_spliceai_out = Channel.value([file('NO_FILE_VCF'), file('NO_FILE_TBI')])
        }

        ch_dbnsfp = Channel.value([
            file(params.dbnsfp),
            file("${params.dbnsfp}.tbi")
            ])

        VEP_GERMLINE_SNV(ch_vcf,
                ch_tbi,
                ch_spliceai_out,
                ch_dbnsfp,
                ch_fasta
                )
        
        FILTER_VEP(VEP_GERMLINE_SNV.out)

}