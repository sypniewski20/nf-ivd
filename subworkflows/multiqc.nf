include { MULTIQC } from '../modules/seqQC.nf'

// ==========================
// MULTIQC
// ==========================
workflow multiqc_workflow {

    take:
        fastqc
        flagstat
        mosdepth
        vcf_stats

    main:
        MULTIQC(
            fastqc.collect(),
            flagstat.collect(),
            mosdepth.collect(),
            vcf_stats.collect()
        )
}