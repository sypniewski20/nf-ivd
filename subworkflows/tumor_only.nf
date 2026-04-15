
include {
    MUTECT2_SOMATIC_ONLY;
    MUTECT2_ORIENTATION_MODEL;
    MUTECT2_ARTIFACT_METRICS;
    MUTECT2_FILTER;
    MUTECT2_EXTRACT_GT;
} from "../modules/Mutect2.nf"


workflow tumor_only_workflow {

    take:
        ch_bam
        params


    main:

        // ============================================================
        // CONSTANTS
        // ============================================================
        def fasta_bundle = tuple(params.fasta, params.fai, params.fasta_dict)

        def snv_bundle   = tuple(params.snv_resource, params.snv_tbi)

        def intervals    = file(params.intervals)

        def pon_bundle = tuple(params.pon, params.pon_tbi)

        MUTECT2_SOMATIC_ONLY(
            ch_bam,
            fasta_bundle,
            snv_bundle,
            pon_bundle,
            params.interval_padding,
            intervals
        )

        // ============================================================
        // 2. ORIENTATION MODEL
        // ============================================================
        def rom = MUTECT2_ORIENTATION_MODEL(
            mutect.f1r2
        )


        // ============================================================
        // 3. ARTIFACT METRICS (independent per sample)
        // ============================================================
        def artifacts = MUTECT2_ARTIFACT_METRICS(
            ch_bam,
            fasta_bundle
        )


        // ============================================================
        // 4. FILTERING (requires ROM)
        // ============================================================
        def filtered = MUTECT2_FILTER(
            mutect.vcf,
            fasta_bundle,
            rom
        )


        // ============================================================
        // 5. EXTRACT GENOTYPES
        // ============================================================
        def gt = MUTECT2_EXTRACT_GT(
            filtered
        )


    emit:
        vcf      = filtered
        gt       = gt
        artifacts = artifacts
}