nextflow.enable.dsl=2

include {
    HAPY_GERMLINE_EVAL
    HAPY_STRATIFIED_EVAL
    HAPY_EXTRACT_METRICS
    AGGREGATE_JSON_REPORT
    GENERATE_CLINICAL_PDF
} from "../modules/HapPy.nf"

workflow GERMLINE_BENCHMARK {

    take:
        ch_vcf
        params

    main:

        def fasta = file(params.fasta)

        // ------------------------------------------------------------
        // GIAB SAMPLE LIST (REGISTRY-DRIVEN)
        // ------------------------------------------------------------
        def giab_samples = ["HG002","HG003","HG004","HG005","HG006","HG007"]


        // ------------------------------------------------------------
        // CORE HAP.PY EVALUATION (GERMLINE)
        // ------------------------------------------------------------
        def evals = giab_samples
            .collect { sample ->

                def truth_vcf = file(params.truth[sample].vcf)
                def truth_bed = file(params.truth[sample].bed)

                HAPY_GERMLINE_EVAL(
                    tuple(sample, ch_vcf),
                    tuple(truth_vcf, truth_bed),
                    fasta
                )
            }


        // ------------------------------------------------------------
        // METRICS EXTRACTION
        // ------------------------------------------------------------
        def metrics = evals
            .map { e -> HAPY_EXTRACT_METRICS(e) }


        // ------------------------------------------------------------
        // STRATIFIED ANALYSIS (HG002 ONLY - GIAB STANDARD)
        // ------------------------------------------------------------
        def hg002_truth = file(params.truth.HG002.vcf)
        def hg002_bed   = file(params.truth.HG002.bed)
        def strat_dir   = file(params.strat_dir)

        def strat_eval = HAPY_STRATIFIED_EVAL(
            tuple("HG002", ch_vcf),
            tuple(hg002_truth, hg002_bed),
            strat_dir,
            fasta
        )


        // ------------------------------------------------------------
        // CLINICAL BUNDLE OUTPUT
        // ------------------------------------------------------------
        def json = AGGREGATE_JSON_REPORT(metrics)

        GENERATE_CLINICAL_PDF(json)
}