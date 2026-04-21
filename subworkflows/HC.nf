// ============================================================
// GERMLINE VARIANT CALLING SUBWORKFLOW (WES + WGS SUPPORT)
// ============================================================

include {
    CALIBRATE_DRAGSTR_MODEL;
    GVCF_HAPLOTYPE_CALLER;
    GENOMICSDB_IMPORT;
    GENOTYPE_GVCF;
    VARIANT_FILTERING;
    CALCULATE_POSTERIORS;
    HAPLOTYPE_CALLER_EXTRACT_GT
} from "../modules/HaplotypeCaller.nf"

workflow hc_workflow {

    take:
        ch_bam       // tuple(val(sample), path(bam), path(bai))

    main:
        
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

        ch_dragstr = CALIBRATE_DRAGSTR_MODEL(
            ch_bam, 
            ch_fasta, 
            file(params.bed), 
            params.interval_padding
         )
        // ============================================================
        // STEP 2: CALLING 
        // ============================================================
        
        ch_gvcfs_out = GVCF_HAPLOTYPE_CALLER(
            ch_bam.join(ch_dragstr),
            ch_fasta,
            file(params.bed),
            params.interval_padding
        )

        ch_input_vcf = ch_gvcfs_out.vcf.collect()
        ch_input_tbi = ch_gvcfs_out.tbi.collect()

        ch_db = GENOMICSDB_IMPORT(
            ch_input_vcf, 
            ch_input_tbi, 
            ch_fasta, 
            file(params.bed))

        ch_raw_vcf = GENOTYPE_GVCF(ch_db, ch_fasta, file(params.bed))

        ch_filtered = VARIANT_FILTERING(ch_raw_vcf, ch_fasta)

        if (params.pedigree) {
            ch_post = CALCULATE_POSTERIORS(
                ch_filtered.vcf, 
                ch_filtered.tbi, 
                file(params.pedigree)
            )

            HAPLOTYPE_CALLER_EXTRACT_GT(ch_post.vcf, ch_post.tbi)
        }

        HAPLOTYPE_CALLER_EXTRACT_GT(ch_filtered.vcf, ch_filtered.tbi)

    emit:
        hc_vcf = ch_filtered.vcf
        stats = ch_filtered.stats
}