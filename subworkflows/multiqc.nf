include { MULTIQC } from '../modules/seqQC.nf'

// ==========================
// MULTIQC
// ==========================
workflow multiqc {

    take:
        fastqc
        flagstat
        mosdepth

    main:
        MULTIQC(
            fastqc.collect(),
            flagstat.collect(),
            mosdepth.collect()
        )
}