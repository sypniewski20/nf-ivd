def Read_samplesheet(samplesheet) {
    // We return the channel object so main.nf can receive it
    return Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            // Use tuple for clarity in DSL2
            tuple(
                row.SM,
                row.ID ?: "${row.RGSM}_${row.LB}",
                row.LB ?: "unknown_lib",
                row.PL ?: "unknown_platform",
                row.PU ?: "unknown_unit",
                file(row.R1, checkIfExists: true), 
                file(row.R2, checkIfExists: true)
                
            )
        }
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