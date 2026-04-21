// ============================================================
// CLINICAL GATK HAPLOTYPECALLER MODULES
// DRAGEN-MODE + DRAGSTR (WES & WGS VERSIONS)
// ============================================================

def GATK_GLOBAL_ARGS = """
--dragen-mode true \
--native-pair-hmm-threads 16 \
--standard-min-confidence-threshold-for-calling 20
"""

// ------------------------------------------------------------
// 1. CALIBRATION (DRAGSTR)
// ------------------------------------------------------------

process CALIBRATE_DRAGSTR_MODEL {
    label 'medium'
    tag "${sample}"
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/HC/dragstr", mode: 'copy'

    input:
        tuple val(sample), path(bam), path(bai)
        tuple path(fasta), path(fai), path(fasta_dict), path(str_table)
        path(bed)
        val(interval_padding)

    output:
        tuple val(sample), path("${sample}_dragstr_model.txt")

    script:
        def padding = (params.seq_type == 'WES') ? params.interval_padding : 0
    """
    gatk CalibrateDragstrModel \
        -R ${fasta} \
        -I ${bam} \
        -L ${bed} \
        --interval-set-rule INTERSECTION \
        --interval-padding ${padding} \
        -str ${str_table} \
        -O ${sample}_dragstr_model.txt
    """
}

// ------------------------------------------------------------
// 2. HAPLOTYPE CALLING (GVCF)
// ------------------------------------------------------------

process GVCF_HAPLOTYPE_CALLER {
    label 'xlarge'
    tag "${sample}"
    label 'gatk'

    input:
        tuple val(sample), path(bam), path(bai), path(dragstr)
        tuple path(fasta), path(fai), path(fasta_dict)
        path(bed)
        val(interval_padding)

    output:
        tuple val(sample), path("${sample}.g.vcf.gz"), emit: vcf
        tuple val(sample), path("${sample}.g.vcf.gz.tbi"), emit: tbi
    script:
        def padding = (params.seq_type == 'WES') ? params.interval_padding : 0
    """
    gatk --java-options "-Xmx32g" HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample}.g.vcf.gz \
        -L ${bed} \
        --interval-set-rule INTERSECTION \
        --interval-padding ${padding} \
        --dragstr-params-path ${dragstr} \
        ${GATK_GLOBAL_ARGS} \
        -ERC GVCF
    """
}

// ------------------------------------------------------------
// 3. JOINT GENOTYPING
// ------------------------------------------------------------

process GENOMICSDB_IMPORT {
    label 'xlarge'
    label 'gatk'

    input:
        path(gvcfs)
        path(tbis)
        tuple path(fasta), path(fai), path(fasta_dict)
        path(bed)

    output:
        path("genomics_db")

    script:
        def input_files = gvcfs.collect { "-V $it" }.join(' ')
    """

    gatk --java-options "-Xmx32g -Xms32g" GenomicsDBImport \
        --genomicsdb-workspace-path genomics_db \
        -R ${fasta} \
        -L ${bed} \
        ${input_files} \
        --tmp-dir . 
    """
}

process GENOTYPE_GVCF {
    label 'large'
    label 'gatk'

    input:
        path(gendb)
        tuple path(fasta), path(fai), path(fasta_dict)
        path(bed)

    output:
        tuple path("HC_joint.vcf.gz"), path("HC_joint.vcf.gz.tbi")
    script:
    """
    gatk --java-options "-Xmx16g" GenotypeGVCFs \
        -R ${fasta} \
        -V gendb://${gendb} \
        -L ${bed} \
        -O HC_joint.vcf.gz

    tabix -p vcf HC_joint.vcf.gz
    """
}

// ------------------------------------------------------------
// 4. GATHER & FILTERING
// ------------------------------------------------------------

process VARIANT_FILTERING {
    label 'medium'
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/HC/filtered", mode: 'copy'

    input:
        tuple path(vcf), path(tbi)
        tuple path(fasta), path(fai), path(fasta_dict)
    output:
        path("HC_filtered_norm.vcf.gz"), emit: vcf
        path("HC_filtered_norm.vcf.gz.tbi"), emit: tbi
        path("HC_filtered_norm.vcf.gz.stats"), emit: stats
        path("HC_filtered_norm.vcf.gz.md5"), emit: md5
    script:
    """
    
    gatk VariantFiltration \
        -R ${fasta} \
        -V ${vcf} \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0" \
        --filter-name "GATK_HARD_FILTER" \
        -O - | \
    bcftools norm -a --atom-overlaps . -m - -f ${fasta} -Ou | \
    bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' -Ou | \
    bcftools +fill-tags -Ou -- -t AF,AC | \
    bcftools sort -Oz -o HC_filtered_norm.vcf.gz

    tabix -p vcf HC_filtered_norm.vcf.gz

    bcftools stat HC_filtered_norm.vcf.gz > HC_filtered_norm.vcf.gz.stats
    md5sum HC_filtered_norm.vcf.gz > HC_filtered_norm.vcf.gz.md5

    """
}

process CALCULATE_POSTERIORS {
    label 'gatk'
    label 'medium'
    publishDir "${params.outfolder}/${params.runID}/HC/posteriors", mode: 'copy'

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
    label 'tiny'
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/HC/genotypes", mode: 'copy'

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
        -F CHROM -F POS -F ID -F REF -F ALT -GF GT -GF AD -GF DP -GF GQ \
        -O ${vcf.baseName}_gt.table

    md5sum ${vcf.baseName}_gt.table > ${vcf.baseName}_gt.table.md5
    """
}