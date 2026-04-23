// ============================================================
// CLINICAL GATK HAPLOTYPECALLER MODULES
// DRAGEN-MODE + DRAGSTR (WES & WGS VERSIONS)
// ============================================================

def GATK_GLOBAL_ARGS = "--dragen-mode true \
--native-pair-hmm-threads 16 \
--standard-min-confidence-threshold-for-calling 20"

def GLOBAL_JAVA_OPTS = "-Xmx32g -XX:+UseParallelGC -XX:ParallelGCThreads=16"

def padding = (params.seq_type == 'WES') ? params.interval_padding : 0


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
        tuple path(fasta), path(fai), path(fasta_dict), path(str_table)
        path(bed)
        val(interval_padding)

    output:
        path("${sample}.g.vcf.gz"), emit: vcf
        path("${sample}.g.vcf.gz.tbi"), emit: tbi
    script:
    """
    gatk --java-options "${GLOBAL_JAVA_OPTS}" HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample}.g.vcf.gz \
        -L ${bed} \
        --interval-set-rule INTERSECTION \
        --interval-padding ${padding} \
        --dragstr-params-path ${dragstr} \
        --smith-waterman FASTEST_AVAILABLE \
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
    tag "${chrom}"
    input:
        val(chrom)
        path(gvcfs)
        path(tbis)
        tuple path(fasta), path(fai), path(fasta_dict), path(str_table)

    output:
        tuple val(chrom), path("genomicsdb_${chrom}")

    script:
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
    label 'large'
    label 'gatk'
    tag "${chrom}"
    input:
        tuple val(chrom), path(gendb)
        tuple path(fasta), path(fai), path(fasta_dict), path(str_table)

    output:
        path("HC_${chrom}.vcf.gz"), emit: vcf
        path("HC_${chrom}.vcf.gz.tbi"), emit: tbi
    script:
    """
    gatk --java-options "${GLOBAL_JAVA_OPTS}" GenotypeGVCFs \
        -R ${fasta} \
        -V gendb://${gendb} \
        -L ${chrom} \
        -O HC_${chrom}.vcf.gz
    """
}

// ------------------------------------------------------------
// 4. GATHER & FILTERING
// ------------------------------------------------------------

process COLLECT_AND_VARIANT_FILTERING {
    label 'medium'
    label 'gatk'
    publishDir "${params.outfolder}/${params.runID}/HC/filtered", mode: 'copy'

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
        -V ${vcf} \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0" \
        --filter-name "GATK_HARD_FILTER" \
        -O HC_raw_tagged.vcf.gz

    bcftools norm -a --atom-overlaps . -m - -f ${fasta} HC_temp_tagged.vcf.gz -Ou | \
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
        -F CHROM -F POS -F ID -F REF -F ALT -GF GT -GF DP \
        -O ${vcf.baseName}_gt.table

    md5sum ${vcf.baseName}_gt.table > ${vcf.baseName}_gt.table.md5
    """
}