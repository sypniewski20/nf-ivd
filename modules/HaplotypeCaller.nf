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

process WES_CALIBRATE_DRAGSTR_MODEL {
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

process WGS_CALIBRATE_DRAGSTR_MODEL {
    label 'medium'
    tag "${sample}"
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/HC/dragstr", mode: 'copy'

    input:
        tuple val(sample), path(bam), path(bai)
        tuple path(fasta), path(fai), path(fasta_dict), path(str_table)

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

// ------------------------------------------------------------
// 2. HAPLOTYPE CALLING (GVCF)
// ------------------------------------------------------------

process WES_GVCF_HAPLOTYPE_CALLER {
    label 'xlarge'
    tag "${sample}:WES"
    label 'gatk'

    input:
        tuple val(sample), path(bam), path(bai), path(dragstr)
        tuple path(fasta), path(fai), path(fasta_dict)
        path(bed)
        val(interval_padding)

    output:
        tuple val(sample), val("exome"), path("${sample}_exome.g.vcf.gz"), emit: vcf
        tuple val(sample), val("exome"), path("${sample}_exome.g.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk --java-options "-Xmx32g" HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample}_exome.g.vcf.gz \
        -L ${bed} \
        --interval-set-rule INTERSECTION \
        --interval-padding ${interval_padding} \
        --dragstr-params-path ${dragstr} \
        ${GATK_GLOBAL_ARGS} \
        -ERC GVCF
    """
}

process WGS_GVCF_HAPLOTYPE_CALLER {
    label 'large'
    tag "${sample}:${interval_name}"
    label 'gatk'

    input:
        tuple val(sample), path(bam), path(bai), path(dragstr), val(interval_name)
        tuple path(fasta), path(fai), path(fasta_dict)

    output:
        tuple val(sample), val(interval_name), path("${sample}_${interval_name}.g.vcf.gz"), emit: vcf
        tuple val(sample), val(interval_name), path("${sample}_${interval_name}.g.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk --java-options "-Xmx16g" HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample}_${interval_name}.g.vcf.gz \
        -L ${interval_name} \
        --dragstr-params-path ${dragstr} \
        ${GATK_GLOBAL_ARGS} \
        -ERC GVCF
    """
}

// ------------------------------------------------------------
// 3. JOINT GENOTYPING
// ------------------------------------------------------------

process WES_GENOMICSDB_IMPORT {
    label 'xlarge'
    tag "WES_DB"
    label 'gatk'

    input:
        tuple path(gvcfs), path(tbis)
        tuple path(fasta), path(fai), path(fasta_dict)
        path(bed)

    output:
        tuple val("exome"), path("wes_genomics_db")

    script:
    """
    rm -rf wes_genomics_db
    V_ARGS=\$(echo ${gvcfs} | tr ' ' '\\n' | sed 's/^/-V /' | tr '\\n' ' ')

    gatk --java-options "-Xmx32g -Xms32g" GenomicsDBImport \
        --genomicsdb-workspace-path wes_genomics_db \
        -R ${fasta} \
        -L ${bed} \
        \$V_ARGS \
        --tmp-dir . \
        --reader-threads 2
    """
}

process WGS_GENOMICSDB_IMPORT {
    label 'xlarge'
    tag "${interval_name}"
    label 'gatk'

    input:
        tuple val(interval_name), path(gvcfs), path(tbis)
        tuple path(fasta), path(fai), path(fasta_dict)

    output:
        tuple val(interval_name), path("${interval_name}_db")

    script:
    """
    rm -rf ${interval_name}_db
    V_ARGS=\$(echo ${gvcfs} | tr ' ' '\\n' | sed 's/^/-V /' | tr '\\n' ' ')

    gatk --java-options "-Xmx16g -Xms16g" GenomicsDBImport \
        --genomicsdb-workspace-path ${interval_name}_db \
        -R ${fasta} \
        -L ${interval_name} \
        \$V_ARGS \
        --tmp-dir . \
        --reader-threads 2
    """
}

process WES_GENOTYPE_GVCF {
    label 'large'
    tag "WES_Genotype"
    label 'gatk'

    input:
        path(gendb)
        tuple path(fasta), path(fai), path(fasta_dict)
        path(bed)

    output:
        path("exome_joint.vcf.gz"), emit: vcf

    script:
    """
    gatk --java-options "-Xmx16g" GenotypeGVCFs \
        -R ${fasta} \
        -V gendb://${gendb} \
        -L ${bed} \
        -O exome_joint.vcf.gz
    """
}

process WGS_GENOTYPE_GVCF {
    label 'large'
    tag "${interval_name}"
    label 'gatk'

    input:
        tuple val(interval_name), path(gendb)
        tuple path(fasta), path(fai), path(fasta_dict)

    output:
        path("${interval_name}.vcf.gz"), emit: vcf

    script:
    """
    gatk --java-options "-Xmx16g" GenotypeGVCFs \
        -R ${fasta} \
        -V gendb://${gendb} \
        -O ${interval_name}.vcf.gz
    """
}

// ------------------------------------------------------------
// 4. GATHER & FILTERING
// ------------------------------------------------------------

process GATHER_VCFS {
    label 'medium'
    tag "${params.runID}"
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/HC/raw", mode: 'copy'

    input:
        path(vcf_list) 

    output:
        path("${params.runID}_raw.vcf.gz"), emit: vcf
        path("${params.runID}_raw.vcf.gz.tbi"), emit: tbi

    script:
    """
    V_INPUTS=\$(echo ${vcf_list} | tr ' ' '\\n' | sed 's/^/-I /' | tr '\\n' ' ')

    gatk GatherVcfs \
        \$V_INPUTS \
        -O ${params.runID}_raw.vcf.gz

    gatk IndexFeatureFile -I ${params.runID}_raw.vcf.gz
    """
}

process VARIANT_FILTERING {
    label 'medium'
    tag "${vcf.baseName}"
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/HC/filtered", mode: 'copy'

    input:
        path(vcf)
        path(tbi)
        tuple path(fasta), path(fai), path(fasta_dict)

    output:
        path("${params.runID}_filtered.vcf.gz"), emit: vcf
        path("${params.runID}_filtered.vcf.gz.tbi"), emit: tbi
        path("${params.runID}_filtered.vcf.gz.stats"), emit: stats
        path("${params.runID}_filtered.vcf.gz.md5"), emit: md5
    script:
    """
    # Applying GATK hard filters as a baseline
    gatk VariantFiltration \
        -R ${fasta} \
        -V ${vcf} \
        -O ${params.runID}_filtered.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0" \
        --filter-name "GATK_HARD_FILTER"

    bcftools stat ${params.runID}_filtered.vcf.gz > ${params.runID}_filtered.vcf.gz.stats
    md5sum ${params.runID}_filtered.vcf.gz > ${params.runID}_filtered.vcf.gz.md5

    """
}

process CALCULATE_POSTERIORS {
    tag "${vcf.baseName}"
    label 'gatk'
    label 'medium'
    publishDir "${params.outfolder}/${params.runID}/HC/posteriors", mode: 'copy'

    input:
        path(vcf)
        path(vcf_index)
        path(pedigree)

    output:
        path("${params.runID}_posteriors.vcf.gz"), emit: vcf
        path("${params.runID}_posteriors.vcf.gz.tbi"), emit: tbi
        path("${params.runID}_posteriors.vcf.gz.md5"), emit: md5
    script:
    """
    gatk CalculateGenotypePosteriors \
        -V ${vcf} \
        -ped ${pedigree} \
        -O ${params.runID}_posteriors.vcf.gz \
        ${GATK_GLOBAL_ARGS}

    tabix -p vcf ${params.runID}_posteriors.vcf.gz
    md5sum ${params.runID}_posteriors.vcf.gz > ${params.runID}_posteriors.vcf.gz.md5

    """
}

process HAPLOTYPE_CALLER_EXTRACT_GT {
    label 'tiny'
    tag "${vcf.baseName}"
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/HC/genotypes", mode: 'copy'

    input:
        path(vcf)

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