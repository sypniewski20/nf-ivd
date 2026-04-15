include { Read_samplesheet } from '../modules/functions.nf'

// ==========================
// PREPROCESSING
// ==========================
workflow preprocessing_workflow {

    take:
        samplesheet

    main:
        Read_samplesheet(samplesheet)

    emit:
        fq = Read_samplesheet.out
}