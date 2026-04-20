// ============================================================
// GATK gCNV subworkflow
// ============================================================
 include { PREPROCESS_INTERVALS; 
           ANNOTATE_INTERVALS; 
           COLLECT_READ_COUNTS; 
           FILTER_INTERVALS;
           DETERMINE_PLOIDY;
           GERMLINE_CNV_CALLER;
           POSTPROCESS_CALLS } from '../modules/gCNV.nf'


workflow gCNV_workflow {
    take:
        ch_bam
    
    main:
        ch_fasta = Channel.value([
            file(params.fasta),
            file("${params.fasta}.fai")
            ])

        bed = Channel.value(file(params.bed))
        ploidy_priors = Channel.value(file(params.ploidy_priors))

        // 1. Preprocess Intervals
        preprocessed_intervals = PREPROCESS_INTERVALS(bed, ch_fasta)
        
        // 2. Annotate Intervals
        annotated_intervals = ANNOTATE_INTERVALS(preprocessed_intervals, ch_fasta)
        
        // 3. Collect Read Counts
        hdf5_counts = COLLECT_READ_COUNTS(ch_bam, preprocessed_intervals, ch_fasta)
        

        cohort_hdf5 = hdf5_count.collect()

        // 4. Filter Intervals
        filtered_intervals = FILTER_INTERVALS(cohort_hdf5, preprocessed_intervals, annotated_intervals)
        
        // 5. Determine Ploidy
        ploidy_calls = DETERMINE_PLOIDY(cohort_hdf5, filtered_intervals, annotated_intervals, ploidy_priors)
        
        // 6. Germline CNV Calling
        cnv_calls = GERMLINE_CNV_CALLER(cohort_hdf5, filtered_intervals, annotated_intervals, ploidy_calls)
        
        // 7. Postprocess Calls (Add Annotations + VCF Formatting)
        // final_calls = POSTPROCESS_CALLS(cnv_calls, annotated_intervals, ch_fasta)
}