process CALIBRATE_DRAGSTR_MODEL {
    tag "${sample}"
    label 'gatk'
    label 'medium'
    publishDir "${params.outfolder}/${params.runID}/HC/dragstr", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(bam), path(bai)
        tuple path(fasta), path(fai), path(fasta_dict), path(str_table)
        path(bed)
        val(interval_padding)

    output:
        tuple val(sample), path("${sample}_dragstr_model.txt")

    script:
        """
        gatk CalibrateDragstrModel \
            -R ${fasta} \
            -I ${bam} \
            -L ${bed} \
            --interval-set-rule INTERSECTION \
            --interval-padding ${interval_padding} \
            -str ${str_table} \
            -O ${sample}_dragstr_model.txt
        """
}

process GVCF_HAPLOTYPE_CALLER {
    tag "${sample}"
    label 'gatk'
    label 'xlarge'

    input:
        tuple val(sample), path(bam), path(bai), path(dragstr_model)
        tuple path(fasta), path(fai), path(fasta_dict), path(str_table)
        path(bed)
        val(interval_padding)

    output:
        path("${sample}.g.vcf.gz"), emit: vcf
        path("${sample}.g.vcf.gz.tbi"), emit: tbi

    script:
        def GATK_GLOBAL_ARGS = "--dragen-mode true --native-pair-hmm-threads 16 --standard-min-confidence-threshold-for-calling 20"
        def GLOBAL_JAVA_OPTS = "-Xmx32g -XX:+UseParallelGC -XX:ParallelGCThreads=16"
        def dragstr_arg = dragstr_model.name != "NO_FILE" ? "--dragstr-params ${dragstr_model}" : ""
        """
        gatk --java-options "${GLOBAL_JAVA_OPTS}" HaplotypeCaller \
            -R ${fasta} \
            -I ${bam} \
            -O ${sample}.g.vcf.gz \
            -L ${bed} \
            --interval-set-rule INTERSECTION \
            --interval-padding ${interval_padding} \
            ${dragstr_arg} \
            --smith-waterman FASTEST_AVAILABLE \
            ${GATK_GLOBAL_ARGS} \
            -ERC GVCF
        """
}

process GENOMICSDB_IMPORT {
    tag "${chrom}"
    label 'gatk'
    label 'xlarge'

    input:
        val(chrom)
        path(gvcfs)
        path(tbis)
        tuple path(fasta), path(fai), path(fasta_dict), path(str_table)

    output:
        tuple val(chrom), path("genomicsdb_${chrom}")

    script:
        def GLOBAL_JAVA_OPTS = "-Xmx32g -XX:+UseParallelGC -XX:ParallelGCThreads=16"
        def input_files = gvcfs.collect { "-V $it" }.join(' ')
        """
        gatk --java-options "${GLOBAL_JAVA_OPTS}" GenomicsDBImport \
            --genomicsdb-workspace-path genomicsdb_${chrom} \
            -R ${fasta} \
            -L ${chrom} \
            ${input_files} \
            --tmp-dir . \
            --batch-size ${gvcfs.size()} \
            --bypass-feature-reader
        """
}

process GENOTYPE_GVCF {
    tag "${chrom}"
    label 'gatk'
    label 'large'

    input:
        tuple val(chrom), path(gendb)
        tuple path(fasta), path(fai), path(fasta_dict), path(str_table)

    output:
        path("HC_${chrom}.vcf.gz"), emit: vcf
        path("HC_${chrom}.vcf.gz.tbi"), emit: tbi

    script:
        def GLOBAL_JAVA_OPTS = "-Xmx32g -XX:+UseParallelGC -XX:ParallelGCThreads=16"
        """
        gatk --java-options "${GLOBAL_JAVA_OPTS}" GenotypeGVCFs \
            -R ${fasta} \
            -V gendb://${gendb} \
            -L ${chrom} \
            -O HC_${chrom}.vcf.gz
        """
}

process COLLECT_AND_VARIANT_FILTERING {
    label 'gatk'
    label 'medium'
    publishDir "${params.outfolder}/${params.runID}/HC/filtered", mode: 'copy', overwrite: true

    input:
        path(vcf)
        path(tbi)
        tuple path(fasta), path(fai), path(fasta_dict), path(str_table)

    output:
        path("HC_filtered_norm.vcf.gz"), emit: vcf
        path("HC_filtered_norm.vcf.gz.tbi"), emit: tbi
        path("HC_filtered_norm.vcf.gz.stats"), emit: stats
        path("HC_filtered_norm.vcf.gz.md5"), emit: md5
        path("HC_raw_tagged.vcf.gz"), emit: raw_vcf
        path("HC_raw_tagged.vcf.gz.tbi"), emit: raw_tbi

    script:
        """
        bcftools concat -a -Oz -o HC_raw.vcf.gz ${vcf}
        tabix -p vcf HC_raw.vcf.gz

        gatk VariantFiltration \
            -R ${fasta} \
            -V HC_raw.vcf.gz \
            --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0" \
            --filter-name "GATK_HARD_FILTER" \
            -O HC_raw_tagged.vcf.gz

        bcftools norm -a --atom-overlaps . -m - -f ${fasta} HC_raw_tagged.vcf.gz -Ou | \
        bcftools view -f PASS -Ou | \
        bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' -Ou | \
        bcftools +fill-tags -Ou -- -t AF,AC | \
        bcftools sort -Oz -o HC_filtered_norm.vcf.gz

        tabix -p vcf HC_filtered_norm.vcf.gz

        bcftools stats HC_filtered_norm.vcf.gz > HC_filtered_norm.vcf.gz.stats
        md5sum HC_filtered_norm.vcf.gz > HC_filtered_norm.vcf.gz.md5
        """
}

process CALCULATE_POSTERIORS {
    label 'gatk'
    label 'medium'
    publishDir "${params.outfolder}/${params.runID}/HC/posteriors", mode: 'copy', overwrite: true

    input:
        path(vcfs)
        path(tbis)
        path(pedigree)

    output:
        path("HC_posteriors.vcf.gz"), emit: vcf
        path("HC_posteriors.vcf.gz.tbi"), emit: tbi
        path("HC_posteriors.vcf.gz.md5"), emit: md5

    script:
        def input_files = vcfs.collect { "-V $it" }.join(' ')
        """
        gatk CalculateGenotypePosteriors \
             ${input_files} \
            -ped ${pedigree} \
            -O HC_posteriors.vcf.gz

        tabix -p vcf HC_posteriors.vcf.gz
        md5sum HC_posteriors.vcf.gz > HC_posteriors.vcf.gz.md5
        """
}

process HAPLOTYPE_CALLER_EXTRACT_GT {
    label 'gatk'
    label 'tiny'
    publishDir "${params.outfolder}/${params.runID}/HC/genotypes", mode: 'copy', overwrite: true

    input:
        path(vcf)
        path(tbi)

    output:
        path("${vcf.baseName}_gt.table")
        path("${vcf.baseName}_gt.table.md5")

    script:
        """
        gatk VariantsToTable \
            -V ${vcf} \
            -F CHROM -F POS -F ID -F REF -F ALT -GF GT -GF DP \
            -O ${vcf.baseName}_gt.table

        md5sum ${vcf.baseName}_gt.table > ${vcf.baseName}_gt.table.md5
        """
}