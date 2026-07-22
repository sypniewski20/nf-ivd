// workflows/cnvkit_workflow.nf

include { CNVKIT_TARGET;
          CNVKIT_ANTITARGET;
          CNVKIT_COVERAGE;
          CNVKIT_REFERENCE;
          CNVKIT_FIX;
          CNVKIT_SEGMENT;
          CNVKIT_CALL;
          CNVKIT_SCATTER;
          CNVKIT_DIAGRAM;
          CNVKIT_TO_VCF;
          MERGE_VCF
        } from '../modules/cnvkit.nf'

workflow cnvkit_workflow {
    take:
        ch_bam

    main:
        ch_fasta = Channel.value([
            file(params.fasta),
            file("${params.fasta}.fai")
        ])

        bed = Channel.value(file(params.bed))

        CNVKIT_TARGET(bed)
        CNVKIT_ANTITARGET(CNVKIT_TARGET.out)

        CNVKIT_COVERAGE(
            ch_bam,
            CNVKIT_TARGET.out,
            CNVKIT_ANTITARGET.out
        )

        CNVKIT_REFERENCE(ch_fasta, CNVKIT_TARGET.out, CNVKIT_ANTITARGET.out)

        CNVKIT_FIX(
            CNVKIT_COVERAGE.out,
            CNVKIT_REFERENCE.out
        )

        CNVKIT_SEGMENT(CNVKIT_FIX.out.cnr)
        CNVKIT_CALL(CNVKIT_SEGMENT.out.cns)

        viz_ch = CNVKIT_FIX.out.cnr
            .join(CNVKIT_CALL.out.cns, by: 0)

        CNVKIT_SCATTER(viz_ch)
        CNVKIT_DIAGRAM(viz_ch)
        CNVKIT_TO_VCF(CNVKIT_CALL.out.cns)

        // Renamed variables to avoid shadowing the closure parameters
        ch_vcf_list = CNVKIT_TO_VCF.out
                            .map{ vcf_file, tbi_file -> vcf_file }
                            .collect()
                            
        ch_tbi_list = CNVKIT_TO_VCF.out
                            .map{ vcf_file, tbi_file -> tbi_file }
                            .collect()

        MERGE_VCF(ch_vcf_list, ch_tbi_list)
}