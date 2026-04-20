process PREPROCESS_INTERVALS {
    label 'medium'
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/gCNV", mode: 'copy'

    input:
        path(bed)
        tuple path(fasta), path(fai)

    output:
        path("preprocessed_intervals.interval_list")

    script:
    """
    gatk PreprocessIntervals \
      -L ${bed} \
      -R ${fasta} \
      --bin-length 0 \
      --padding ${params.interval_padding} \
      -O preprocessed_intervals.interval_list
    """
}

process ANNOTATE_INTERVALS {
    label 'medium'
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/gCNV", mode: 'copy'

    input:
        path(preprocessed_intervals)
        tuple path(fasta), path(fai) 

    output:
        path("annotated_intervals.tsv")
    script:
    """

    gatk AnnotateIntervals \
      -L ${preprocessed_intervals} \
      -R ${fasta} \
      -O annotated_intervals.tsv

    """
}

process COLLECT_READ_COUNTS {
    label 'medium'
    label 'gatk'
    tag "${sample}"
    input:
        tuple val(sample), path(bam), path(bai)
        path(preprocessed_intervals)
        tuple path(fasta), path(fai) 

    output:
        path("${sample}.counts.hdf5")
    script:
    """

    gatk CollectReadCounts \
      -I ${bam} \
      -L ${preprocessed_intervals} \
      -R ${fasta} \
      --format HDF5 \
      -O ${sample}.counts.hdf5

    """
}

process FILTER_INTERVALS {
    label 'medium'
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/gCNV", mode: 'copy'

    input:
        path(hdf5_counts)
        path(preprocessed_intervals)
        path(annotated_intervals)
    output:
        path("filtered.interval_list")
    script:
        def input_files = hdf5_counts.collect { "-I $it" }.join(' ')
    
    """
        gatk FilterIntervals \
            -L ${preprocessed_intervals} \
            --annotated-intervals ${annotated_intervals} \
            $input_files \
            -O filtered.interval_list
    """
}

process DETERMINE_PLOIDY {
    publishDir "${params.outfolder}/${params.runID}/gCNV", mode: 'copy'
    label 'medium'
    label 'gatk'
    input:
    path(hdf5_counts)
    path(filtered_intervals)
    path(annotated_intervals)
    path(ploidy_priors)

    output:
    path "ploidy-calls", emit: calls

    script:
    def sorted_files = hdf5_counts.sort { it.name }
    def input_files = sorted_files.collect { "-I $it" }.join(' ')

    """
    gatk DetermineGermlineContigPloidy \
        -L ${filtered_intervals} \
        --annotated-intervals ${annotated_intervals} \
        $input_files \
        --contig-ploidy-priors ${ploidy_priors} \
        --output . \
        --output-prefix ploidy
    """
}

process GERMLINE_CNV_CALLER {
    label 'medium'
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/gCNV", mode: 'copy'

    input:
        path(hdf5_counts)
        path(filtered_intervals)
        path(annotated_intervals)
        path(ploidy_calls)
    output:
        path "cohort_model_output/cohort_run-model", emit: model
        path "cohort_model_output/cohort_run-calls", emit: calls
    script:
        def sorted_files = hdf5_counts.sort { it.name }
        def input_files = sorted_files.collect { "-I $it" }.join(' ')
    
    """
        gatk GermlineCNVCaller \
          --run-mode COHORT \
          -L ${filtered_intervals} \
          --annotated-intervals ${annotated_intervals} \
          $input_files \
          --contig-ploidy-calls ${ploidy_calls} \
          --output cohort_model_output \
          --output-prefix cohort_run
    """
}

process POSTPROCESS_CALLS {
    tag "${sample}"
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/gCNV", mode: 'copy'
    
    input:
    tuple val(index), val(sample)
    path (model_dir)
    path (calls_dir)
    path (ploidy_dir)

    output:
    path "${sample}.intervals.vcf.gz", emit: intervals_vcf
    path "${sample}.segments.vcf.gz",  emit: segments_vcf

    script:
    // We point to the folder AND the prefix 'cohort_run' used in the Caller step
    """
    gatk PostprocessGermlineCNVCalls \
      --model-shard-path ${model_dir} \
      --calls-shard-path ${calls_dir} \
      --all-ploidy-calls ${ploidy_dir} \
      --sample-index ${index} \
      --output-genotyped-intervals ${sample}.intervals.vcf.gz \
      --output-genotyped-segments ${sample}.segments.vcf.gz \
      --output-denoised-copy-ratios ${sample}.denoisedCR.tsv
    """
}