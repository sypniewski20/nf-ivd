def Read_samplesheet(samplesheet) {

    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            [
                row.sampleID,
                file(row.R1, checkIfExists: true),
                file(row.R2, checkIfExists: true),
                [lb: row.LB ?: null]
            ]
        }
        .set { ch_fq }
}

def Read_bam(bam_sheet) {

    Channel
        .fromPath(bam_sheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            [
                row.sampleID,
                file(row.bam, checkIfExists: true),
                file(row.bai, checkIfExists: true)
            ]
        }
        .set { ch_bam }
}