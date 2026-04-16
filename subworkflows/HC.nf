include {
    WES_CALIBRATE_DRAGSTR_MODEL;
    WGS_CALIBRATE_DRAGSTR_MODEL

    WES_GVCF_HAPLOTYPE_CALLER;
    WGS_GVCF_HAPLOTYPE_CALLER

    WES_HAPLOTYPE_CALLER;
    WGS_HAPLOTYPE_CALLER

    WES_GENOTYPE_GVCF;
    WGS_GENOTYPE_GVCF

    WES_GENOMICSDB_IMPORT;
    WGS_GENOMICSDB_IMPORT

    HAPLOTYPE_CALLER_EXTRACT_GT;
    GATHER_VCF;
    VARIANT_FILTERING
} from "../modules/HaplotypeCaller.nf"


def DRAGEN_PROFILE = [
    "--dragen-mode true",
    "--native-pair-hmm-threads 1",
    "--standard-min-confidence-threshold-for-calling 20"
].join(" ")


workflow hc_workflow {

    take:
        ch_bam
        params

    main:

        def isWES = params.seq_type == 'WES'


        // ============================================================
        // 1. DRAGSTR CALIBRATION (MANDATORY FIRST STEP)
        // ============================================================

        def dragstr_model = isWES
            ? WES_CALIBRATE_DRAGSTR_MODEL(
                ch_bam,
                params.fasta,
                params.str_table,
                params.bed,
                params.interval_padding
              )
            : WGS_CALIBRATE_DRAGSTR_MODEL(
                ch_bam,
                params.fasta,
                params.str_table
              )


        // ============================================================
        // 2. GVCF CALLING (DRAGSTR DEPENDENT)
        // ============================================================

        def gvcf = isWES
            ? WES_GVCF_HAPLOTYPE_CALLER(
                ch_bam,
                dragstr_model,
                params.fasta,
                params.interval_padding,
                DRAGEN_PROFILE
              )
            : WGS_GVCF_HAPLOTYPE_CALLER(
                ch_bam,
                dragstr_model,
                params.fasta,
                DRAGEN_PROFILE
              )


        // ============================================================
        // 3. GENOMICSDB
        // ============================================================

        def db = isWES
            ? WES_GENOMICSDB_IMPORT(gvcf.vcf, gvcf.tbi, params.fasta, params.bed)
            : WGS_GENOMICSDB_IMPORT(gvcf.vcf, gvcf.tbi, params.fasta)


        // ============================================================
        // 4. GENOTYPING
        // ============================================================

        def genotype = isWES
            ? WES_GENOTYPE_GVCF(db, params.fasta, params.bed)
            : WGS_GENOTYPE_GVCF(db, params.fasta)


        // ============================================================
        // 5. DIRECT HC (CONSISTENCY LAYER)
        // ============================================================

        def hc = isWES
            ? WES_HAPLOTYPE_CALLER(
                ch_bam,
                dragstr_model,
                params.fasta,
                params.interval_padding,
                DRAGEN_PROFILE
              )
            : WGS_HAPLOTYPE_CALLER(
                ch_bam,
                dragstr_model,
                params.fasta,
                DRAGEN_PROFILE
              )


        // ============================================================
        // 6. MERGE
        // ============================================================

        def merged = genotype.mix(hc.vcf).collect()


        // ============================================================
        // 7. POST PROCESSING
        // ============================================================

        def gathered = GATHER_VCF(merged)

        def filtered = VARIANT_FILTERING(gathered, params.fasta)

        HAPLOTYPE_CALLER_EXTRACT_GT(filtered)
    

    emit:
        hc_vcf = filtered
}