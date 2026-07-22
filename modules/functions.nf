def readSamplesheet(samplesheet) {
    return Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, strip: true)
        .map { row ->
            // --- validate required fields ---
            if (!row.SM) error "Missing SM in samplesheet row: ${row}"
            if (!row.R1) error "Missing R1 path in samplesheet row: ${row}"
            if (!row.R2) error "Missing R2 path in samplesheet row: ${row}"

            def sm = row.SM
            def lb = row.LB ?: "lib_${sm}"
            def id = row.ID ?: "${sm}_${lb}"    
            def pl = row.PL ?: "ILLUMINA"       
            def pu = row.PU ?: "${id}.unknown"  

            def r1 = row.R1 ? file(row.R1, checkIfExists: true) : error("R1 path is null or empty in row: ${row}")
            def r2 = row.R2 ? file(row.R2, checkIfExists: true) : error("R2 path is null or empty in row: ${row}")

            tuple(sm, id, lb, pl, pu, r1, r2)
        }
}

def readBam(bam_sheet) {
    return Channel
        .fromPath(bam_sheet)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            if (!row.sampleID) error "Missing sampleID in BAM samplesheet row: ${row}"
            if (!row.bam)      error "Missing BAM path in samplesheet row: ${row}"
            if (!row.bai)      error "Missing BAI path in samplesheet row: ${row}"

            def bam_file = row.bam ? file(row.bam, checkIfExists: true) : error("BAM path is null or empty in row: ${row}")
            def bai_file = row.bai ? file(row.bai, checkIfExists: true) : error("BAI path is null or empty in row: ${row}")

            tuple(row.sampleID, bam_file, bai_file)
        }
}