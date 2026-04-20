nextflow.enable.dsl=2

include {
    HAPY_GERMLINE_EVAL
    HAPY_STRATIFIED_EVAL
} from "../modules/HapPy.nf"

workflow germline_calibration_workflow {

    take:
        ch_vcf      // tuple(val(runID), path(vcf), path(tbi))

    main:
        // 1. Create a Reference Channel
        ch_fasta = Channel.value([
            file(params.fasta),
            file("${params.fasta}.fai")
        ])

        // 2. Create a Channel from your GIAB Registry
        // This replaces the .collect loop and makes it "Nextflow-native"
        def giab_list = ["HG002","HG003","HG004"]
        
        ch_giab_samples = Channel.fromList(giab_list)
            .map { sample_id ->
                def truth_vcf = file(params.truth[sample_id].vcf)
                def truth_bed = file(params.truth[sample_id].bed)
                return tuple(sample_id, truth_vcf, truth_bed)
            }

        // 3. CORE HAP.PY EVALUATION
        // We combine the input VCF with every GIAB sample truth set
        ch_eval_input = ch_vcf.combine(ch_giab_samples)
        
        // Input: [runID, vcf, tbi, sample_id, truth_vcf, truth_bed]
        ch_eval_results = HAPY_GERMLINE_EVAL(
            ch_eval_input,
            ch_fasta
        )

// ------------------------------------------------------------
        // STRATIFIED ANALYSIS (HG002 ONLY - PARALLELIZED)
        // ------------------------------------------------------------
        def hg002_truth = file(params.truth.HG002.vcf)
        def hg002_bed   = file(params.truth.HG002.bed)
        
        // 1. Create a channel of all BED files in the stratification directory
        ch_strat_beds = Channel.fromPath("${params.strat_dir}/*.bed*")

        // 2. Filter the evaluation results for HG002 and combine with the beds
        // This creates a unique task for every BED file
        ch_strat_input = ch_eval_results.combined_vcf
            .filter { sample, vcf -> sample == "HG002" }
            .combine(ch_strat_beds) 
            .map { sample, vcf, tbi, strat_bed -> 
                tuple(sample, vcf, tbi, hg002_truth, hg002_bed, strat_bed) 
            }

        // 3. Call the updated module
        ch_strat_results = HAPY_STRATIFIED_EVAL(
            ch_strat_input,
            ch_fasta
        )
}