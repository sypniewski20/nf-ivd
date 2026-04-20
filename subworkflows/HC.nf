// ============================================================
// GERMLINE VARIANT CALLING SUBWORKFLOW (WES + WGS SUPPORT)
// ============================================================

include {
    WES_CALIBRATE_DRAGSTR_MODEL;
    WGS_CALIBRATE_DRAGSTR_MODEL;
    WES_GVCF_HAPLOTYPE_CALLER;
    WGS_GVCF_HAPLOTYPE_CALLER;
    WES_GENOMICSDB_IMPORT;
    WGS_GENOMICSDB_IMPORT;
    WES_GENOTYPE_GVCF;
    WGS_GENOTYPE_GVCF;
    GATHER_VCFS;
    VARIANT_FILTERING;
    CALCULATE_POSTERIORS;
    HAPLOTYPE_CALLER_EXTRACT_GT
} from "../modules/HaplotypeCaller.nf"

workflow hc_workflow {

    take:
        ch_bam       // tuple(val(sample), path(bam), path(bai))

    main:
        def isWES = (params.seq_type == 'WES')
        
        // 1. Prepare Reference Channels

        ch_fasta= Channel.value([
            file(params.fasta),
            file("${params.fasta}.fai"),
            file(params.fasta.replace(".fasta", ".dict").replace(".fa", ".dict")),
            file(params.fasta.replace(".fasta", ".str").replace(".fa", ".str"))

        ])

        // ============================================================
        // STEP 1: DRAGSTR CALIBRATION
        // ============================================================
        if (isWES) {
            ch_dragstr = WES_CALIBRATE_DRAGSTR_MODEL(
                ch_bam, 
                ch_fasta, 
                file(params.bed), 
                params.interval_padding
            )
        } else {
            ch_dragstr = WGS_CALIBRATE_DRAGSTR_MODEL(
                ch_bam, 
                ch_fasta
            )
        }

        // ============================================================
        // STEP 2: CALLING (BRANCHING WES vs WGS)
        // ============================================================
        
        if (isWES) {
            // WES PATH: No Sharding (One job per sample)
            ch_wes_input = ch_bam.join(ch_dragstr)
            
            ch_gvcfs_out = WES_GVCF_HAPLOTYPE_CALLER(
                ch_wes_input,
                ch_fasta,
                file(params.bed),
                params.interval_padding
            )

            // Prepare for Joint Genotyping
            ch_db_input = ch_gvcfs_out.vcf
                .map { sample, interval, vcf -> tuple(interval, vcf) }
                .groupTuple()
                .join(
                    ch_gvcfs_out.tbi
                        .map { sample, interval, tbi -> tuple(interval, tbi) }
                        .groupTuple()
                )

            ch_db = WES_GENOMICSDB_IMPORT(ch_db_input, ch_fasta, file(params.bed))
            ch_raw_vcf = WES_GENOTYPE_GVCF(ch_db, ch_fasta, file(params.bed))
            
            // For WES, "gathered" is simply the output of joint genotyping
            ch_gathered = [vcf: ch_raw_vcf, tbi: ch_raw_vcf.map{ it.replace(".vcf.gz", ".vcf.gz.tbi") }]

        } else {
            // WGS PATH: Sharding enabled
            ch_intervals = Channel.fromList(params.intervals_list)
            ch_wgs_input = ch_bam.join(ch_dragstr).combine(ch_intervals)

            ch_gvcfs_out = WGS_GVCF_HAPLOTYPE_CALLER(ch_wgs_input, ch_fasta)

            ch_db_input = ch_gvcfs_out.vcf
                .map { sample, interval, vcf -> tuple(interval, vcf) }
                .groupTuple()
                .join(
                    ch_gvcfs_out.tbi
                        .map { sample, interval, tbi -> tuple(interval, tbi) }
                        .groupTuple()
                )

            ch_db = WGS_GENOMICSDB_IMPORT(ch_db_input, ch_fasta)
            ch_shards = WGS_GENOTYPE_GVCF(ch_db, ch_fasta)

            // Gather all WGS shards into one VCF
            ch_gathered = GATHER_VCFS(ch_shards.collect())
        }

        // ============================================================
        // STEP 3: FILTERING & REFINEMENT (Shared Logic)
        // ============================================================

        ch_filtered = VARIANT_FILTERING(ch_gathered.vcf, ch_gathered.tbi, ch_fasta)

        if (params.pedigree) {
            ch_final = CALCULATE_POSTERIORS(
                ch_filtered.vcf, 
                ch_filtered.tbi, 
                file(params.pedigree)
            )
        } else {
            ch_final = ch_filtered
        }

        HAPLOTYPE_CALLER_EXTRACT_GT(ch_final.vcf)

    emit:
        hc_vcf = ch_final.vcf
        stats = ch_final.stats
}