process PREPROCESS_INTERVALS {
    label 'medium'
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/gCNV", mode: 'copy'

    input:
        path(bed)
        tuple path(fasta), path(fai), path(dict)

    output:
        path("preprocessed_intervals.interval_list")

    script:
        def padding = (params.seq_type == 'WES') ? params.interval_padding : 0
        """
        gatk PreprocessIntervals \
          -L ${bed} \
          -R ${fasta} \
          --bin-length 0 \
          --padding ${padding} \
          -O preprocessed_intervals.interval_list \
          --interval-merging-rule OVERLAPPING_ONLY
        """
}

process ANNOTATE_INTERVALS {
    label 'medium'
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/gCNV", mode: 'copy'

    input:
        path(preprocessed_intervals)
        tuple path(fasta), path(fai), path(dict)

    output:
        path("annotated_intervals.tsv")

    script:
        """
        gatk AnnotateIntervals \
          -L ${preprocessed_intervals} \
          --interval-merging-rule OVERLAPPING_ONLY \
          -R ${fasta} \
          -O annotated_intervals.tsv
        """
}

process COLLECT_READ_COUNTS {
    tag "${sample}"
    label 'medium'
    label 'gatk'

    input:
        tuple val(sample), path(bam), path(bai)
        path(preprocessed_intervals)
        tuple path(fasta), path(fai), path(dict)

    output:
        path("${sample}.counts.hdf5")

    script:
        """
        gatk CollectReadCounts \
          -I ${bam} \
          -L ${preprocessed_intervals} \
          --interval-merging-rule OVERLAPPING_ONLY \
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
            --interval-merging-rule OVERLAPPING_ONLY \
            ${input_files} \
            -O filtered.interval_list
        """
}

process DETERMINE_PLOIDY {
    label 'medium'
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/gCNV", mode: 'copy'

    input:
        path(hdf5_counts)
        path(filtered_intervals)
        path(annotated_intervals)
        path(ploidy_priors)

    output:
        path("ploidy-calls"), emit: calls

    script:
        def sorted_files = hdf5_counts.sort { it.name }
        def input_files = sorted_files.collect { "-I $it" }.join(' ')
        """
        gatk DetermineGermlineContigPloidy \
            -L ${filtered_intervals} \
            --interval-merging-rule OVERLAPPING_ONLY \
            $input_files \
            --contig-ploidy-priors ${ploidy_priors} \
            --output . \
            --output-prefix ploidy
        """
}

process SCATTER_INTERVALS {
    label 'gatk'

    input:
        path(filtered_intervals)

    output:
        path("scatter/*"), emit: scattered_intervals

    script:
        def num_intervals_per_scatter = (params.seq_type == 'WES') ? 15000 : 30000
        """
        mkdir -p scatter
        gatk IntervalListTools \
            --INPUT ${filtered_intervals} \
            --SUBDIVISION_MODE INTERVAL_COUNT \
            --SCATTER_CONTENT ${num_intervals_per_scatter} \
            --OUTPUT scatter
        """
}

process GERMLINE_CNV_CALLER {
    tag "${scattered_intervals.name.tokenize('.').first()}"
    label 'medium'
    label 'gatk'

    input:
        path(hdf5_counts)
        path(scattered_intervals)
        path(annotated_intervals)
        path(ploidy_calls)

    output:
        path("scatter-*-model"), emit: model
        path("scatter-*-calls"), emit: calls

    script:
        def scatter_index = scattered_intervals.name.tokenize('.').first()
        def sorted_files  = hdf5_counts instanceof List ? hdf5_counts.sort { it.name } : [hdf5_counts]
        def input_files   = sorted_files.collect { "-I $it" }.join(' ')
        """
        gatk GermlineCNVCaller \
          --run-mode COHORT \
          -L ${scattered_intervals}/scattered.interval_list \
          --interval-merging-rule OVERLAPPING_ONLY \
          --annotated-intervals ${annotated_intervals} \
          ${input_files} \
          --contig-ploidy-calls ${ploidy_calls} \
          --output . \
          --output-prefix scatter-${scatter_index}
        """
}

process POSTPROCESS_CALLS {
    tag "${sample}"
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/SV/gCNV", mode: 'copy'

    input:
        tuple val(sample), val(index)
        path(model_dir)
        path(calls_dir)
        path(ploidy_dir)

    output:
        tuple path("${sample}.intervals.vcf.gz"), path("${sample}.segments.vcf.gz"), path("${sample}.denoisedCR.tsv")

    script:
        def model_shards = model_dir.collect { "--model-shard-path $it" }.join(' ')
        def call_shards = calls_dir.collect { "--calls-shard-path $it" }.join(' ')
        """
        gatk PostprocessGermlineCNVCalls \
          ${model_shards} \
          ${call_shards} \
          --contig-ploidy-calls ${ploidy_dir} \
          --sample-index ${index} \
          --output-genotyped-intervals ${sample}.intervals.vcf.gz \
          --output-genotyped-segments ${sample}.segments.vcf.gz \
          --output-denoised-copy-ratios ${sample}.denoisedCR.tsv
        """
}