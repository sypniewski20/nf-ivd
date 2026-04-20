
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

    main:

        // ============================================================
        // CONSTANTS
        // ============================================================
        ch_fasta = Channel.value([
        file(params.fasta),
        file("${params.fasta}.fai")
        ])

        snv_bundle   = tuple(params.snv_resource, params.snv_tbi)
        intervals    = file(params.intervals)
        pon_bundle = tuple(params.pon, params.pon_tbi)

        MUTECT2_SOMATIC_ONLY(
            ch_bam,
            ch_fasta,
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
            ch_fasta
        )


        // ============================================================
        // 4. FILTERING (requires ROM)
        // ============================================================
        def filtered = MUTECT2_FILTER(
            mutect.vcf,
            ch_fasta,
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