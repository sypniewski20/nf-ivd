// ============================================================
// CLINICAL GATK HAPLOTYPECALLER WORKFLOW
// DRAGEN-MODE + DRAGSTR (FULLY ENFORCED)
// ============================================================


// ------------------------------------------------------------
// GLOBAL CONSTANTS (APPLIED TO ALL HAPLOTYPECALLER RUNS)
// ------------------------------------------------------------
def GATK_GLOBAL_ARGS = """
--dragen-mode true \
--dragstr-params-path \${dragstr} \
--native-pair-hmm-threads 1 \
--standard-min-confidence-threshold-for-calling 20
"""

// ============================================================
// DRAGSTR CALIBRATION (WES ONLY)
// ============================================================

process WES_CALIBRATE_DRAGSTR_MODEL {
    label 'medium'
    tag "${sample}"
    label 'gatk'

    input:
        tuple val(sample), path(bam), path(bai)
        tuple path(fasta), path(fai), path(fasta_dict)
        path(str_table)
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


// ============================================================
// WES GVCF CALLING
// ============================================================

process WES_GVCF_HAPLOTYPE_CALLER {
    label 'large'
    tag "${sample}:${interval_name}"
    label 'gatk'

    input:
        tuple val(sample), path(bam), path(bai), path(dragstr)
        val(interval_name)
        path(interval)
        tuple path(fasta), path(fai), path(fasta_dict)
        val(interval_padding)

    output:
        tuple val(sample), val(interval_name),
              path("${sample}_${interval_name}.g.vcf.gz"), emit: vcf
        tuple val(sample), val(interval_name),
              path("${sample}_${interval_name}.g.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk --java-options "-Xmx16g" HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample}_${interval_name}.g.vcf.gz \
        -L ${interval} \
        --interval-padding ${interval_padding} \
        ${GATK_GLOBAL_ARGS} \
        -ERC GVCF
    """
}


// ============================================================
// WES DIRECT VCF CALLING
// ============================================================

process WES_HAPLOTYPE_CALLER {
    publishDir "${params.outfolder}/${params.runID}/SNV/haplotype_caller/WES_${sample}", mode: 'copy', overwrite: true
    tag "${sample}:${interval_name}"
    label 'gatk'
    label 'large'

    input:
        tuple val(sample), path(bam), path(bai), path(dragstr)
        val(interval_name)
        path(interval)
        tuple path(fasta), path(fai), path(fasta_dict)
        val(interval_padding)

    output:
        tuple val(sample), val(interval_name),
              path("${sample}_${interval_name}.vcf.gz"), emit: vcf
        tuple val(sample), val(interval_name),
              path("${sample}_${interval_name}.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk --java-options "-Xmx16g" HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample}_${interval_name}.vcf.gz \
        -L ${interval} \
        --interval-padding ${interval_padding} \
        ${GATK_GLOBAL_ARGS}
    """
}


// ============================================================
// WGS DRAGSTR CALIBRATION
// ============================================================

process WGS_CALIBRATE_DRAGSTR_MODEL {
    tag "${sample}"
    label 'gatk'
    label 'medium'

    input:
        tuple val(sample), path(bam), path(bai)
        tuple path(fasta), path(fai), path(fasta_dict)
        path(str_table)

    output:
        tuple val(sample), path("${sample}_dragstr_model.txt")

    script:
    """
    gatk CalibrateDragstrModel \
        -R ${fasta} \
        -I ${bam} \
        -str ${str_table} \
        -O ${sample}_dragstr_model.txt
    """
}


// ============================================================
// WGS GVCF CALLING
// ============================================================

process WGS_GVCF_HAPLOTYPE_CALLER {
    label 'large'
    tag "${sample}:${interval_name}"
    label 'gatk'

    input:
        tuple val(sample), path(bam), path(bai), path(dragstr)
        val(interval_name)
        tuple path(fasta), path(fai), path(fasta_dict)

    output:
        tuple val(sample), val(interval_name),
              path("${sample}_${interval_name}.g.vcf.gz"), emit: vcf
        tuple val(sample), val(interval_name),
              path("${sample}_${interval_name}.g.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk --java-options "-Xmx16g" HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample}_${interval_name}.g.vcf.gz \
        -L ${interval_name} \
        ${GATK_GLOBAL_ARGS} \
        -ERC GVCF
    """
}


// ============================================================
// WGS DIRECT VCF CALLING
// ============================================================

process WGS_HAPLOTYPE_CALLER {
    publishDir "${params.outfolder}/${params.runID}/SNV/haplotype_caller/WGS_${sample}", mode: 'copy', overwrite: true
    tag "${sample}:${interval_name}"
    label 'gatk'
    label 'large'

    input:
        tuple val(sample), path(bam), path(bai), path(dragstr)
        val(interval_name)
        tuple path(fasta), path(fai), path(fasta_dict)

    output:
        tuple val(sample), val(interval_name),
              path("${sample}_${interval_name}.vcf.gz"), emit: vcf
        tuple val(sample), val(interval_name),
              path("${sample}_${interval_name}.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk --java-options "-Xmx16g" HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample}_${interval_name}.vcf.gz \
        -L ${interval_name} \
        ${GATK_GLOBAL_ARGS}
    """
}


// ============================================================
// GENOMICSDB IMPORT
// ============================================================

process WES_GENOMICSDB_IMPORT {
    label 'xlarge'
    tag "${interval_name}"
    label 'gatk'

    input:
        tuple val(interval_name), path(gvcfs)
        tuple val(interval_name), path(tbis)
        tuple path(fasta), path(fai), path(fasta_dict)
        path(intervals)

    output:
        tuple val(interval_name), path("${interval_name}_db")

    script:
    """
    gatk --java-options "-Xmx16g -Xms16g" GenomicsDBImport \
        --genomicsdb-workspace-path ${interval_name}_db \
        -R ${fasta} \
        -L ${intervals} \
        ${gvcfs.collect { "-V ${it}" }.join(' ')}
    """
}


// ============================================================
// GENOTYPING
// ============================================================

process WES_GENOTYPE_GVCF {
    label 'large'
    tag "${interval_name}"
    label 'gatk'

    input:
        tuple val(interval_name), path(gendb)
        tuple path(fasta), path(fai), path(fasta_dict)
        path(intervals)

    output:
        path("${interval_name}.vcf.gz")

    script:
    """
    gatk --java-options "-Xmx16g" GenotypeGVCFs \
        -R ${fasta} \
        -V gendb://${gendb} \
        -O ${interval_name}.vcf.gz \
        -L ${intervals} \
        --merge-input-intervals
    """
}


// ============================================================
// POST PROCESSING
// ============================================================

process GATHER_VCF {
    label 'small'
    tag "gather_vcf"
    label 'gatk'

    publishDir "${params.outfolder}/${params.runID}/SNV/haplotype_caller", mode: 'copy', overwrite: true

    input:
        path(vcf)

    output:
        tuple path("cohort_gather.vcf.gz"), path("cohort_gather.vcf.gz.tbi")

    script:
    """
    bcftools concat ${vcf} -Ou | \
    bcftools sort -Oz -o cohort_gather.vcf.gz

    tabix -p vcf cohort_gather.vcf.gz
    """
}


// ============================================================
// VARIANT FILTERING (CLINICAL NORMALIZATION)
// ============================================================

process VARIANT_FILTERING {
    label 'small'
    tag "variant_filtering"
    label 'gatk'

    publishDir "${params.outfolder}/${params.runID}/SNV/haplotype_caller", mode: 'copy', overwrite: true

    input:
        tuple path(vcf), path(tbi)
        path(fasta)

    output:
        tuple path("cohort_hc_norm.vcf.gz"), path("cohort_hc_norm.vcf.gz.tbi")

    script:
    """
    gatk VariantFiltration \
        -V ${vcf} \
        --filter-expression "QUAL < 10.4139" \
        --filter-name "DRAGENHardQUAL" \
        -O filtered.vcf.gz

    bcftools view -f PASS filtered.vcf.gz -Ou | \
    bcftools norm -a -m -any -f ${fasta} -Oz -o cohort_hc_norm.vcf.gz

    tabix -p vcf cohort_hc_norm.vcf.gz
    """
}


// ============================================================
// GT EXTRACTION
// ============================================================

process HAPLOTYPE_CALLER_EXTRACT_GT {
    label 'small'
    tag "extract_gt"
    label 'gatk'

    publishDir "${params.outfolder}/${params.runID}/SNV/haplotype_caller", mode: 'copy', overwrite: true

    input:
        tuple path(vcf), path(tbi)

    output:
        path("cohort_hc_norm_gt.tsv.gz")

    script:
    """
    echo -e "ID\\tSAMPLE\\tDP\\tAF\\tQUAL\\tGT" | bgzip -c > cohort_hc_norm_gt.tsv.gz

    bcftools query -f "[%ID\\t%SAMPLE\\t%INFO/DP\\t%AF\\t%QUAL\\t%GT\\n]" ${vcf} | \
    bgzip -c >> cohort_hc_norm_gt.tsv.gz
    """
}