// ============================================================
// GERMline HaplotypeCaller Module (CLEAN VERSION)
// ============================================================

// ----------------------------
// GVCF CALLERS
// ----------------------------

process WES_GVCF_HAPLOTYPE_CALLER {

    tag "${sample}:${interval_name}"
    label 'gatk'

    input:
        tuple val(sample), path(bam), path(bai)
        val(interval_name)
        path(interval)
        tuple path(fasta), path(fai), path(fasta_dict)
        val(interval_padding)

    output:
        tuple val(sample), val(interval_name), path("${sample}_${interval_name}.g.vcf.gz"), emit: vcf
        tuple val(sample), val(interval_name), path("${sample}_${interval_name}.g.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample}_${interval_name}.g.vcf.gz \
        -L ${interval} \
        --interval-padding ${interval_padding} \
        -ERC GVCF
    """
}


process WGS_GVCF_HAPLOTYPE_CALLER {

    tag "${sample}:${interval_name}"
    label 'gatk'

    input:
        tuple val(sample), path(bam), path(bai)
        val(interval_name)
        path(interval)
        tuple path(fasta), path(fai), path(fasta_dict)

    output:
        tuple val(sample), val(interval_name), path("${sample}_${interval_name}.g.vcf.gz"), emit: vcf
        tuple val(sample), val(interval_name), path("${sample}_${interval_name}.g.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample}_${interval_name}.g.vcf.gz \
        -L ${interval} \
        -ERC GVCF
    """
}


// ----------------------------
// DIRECT HAPLOTYPE CALLERS
// ----------------------------

process WES_HAPLOTYPE_CALLER {

    tag "${sample}:${interval_name}"
    label 'gatk'

    input:
        tuple val(sample), path(bam), path(bai)
        val(interval_name)
        path(interval)
        tuple path(fasta), path(fai), path(fasta_dict)
        val(interval_padding)

    output:
        tuple val(sample), val(interval_name), path("${sample}_${interval_name}.vcf.gz"), emit: vcf
        tuple val(sample), val(interval_name), path("${sample}_${interval_name}.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample}_${interval_name}.vcf.gz \
        -L ${interval} \
        --interval-padding ${interval_padding}
    """
}


process WGS_HAPLOTYPE_CALLER {

    tag "${sample}:${interval_name}"
    label 'gatk'

    input:
        tuple val(sample), path(bam), path(bai)
        val(interval_name)
        path(interval)
        tuple path(fasta), path(fai), path(fasta_dict)

    output:
        tuple val(sample), val(interval_name), path("${sample}_${interval_name}.vcf.gz"), emit: vcf
        tuple val(sample), val(interval_name), path("${sample}_${interval_name}.vcf.gz.tbi"), emit: tbi

    script:
    """
    gatk HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample}_${interval_name}.vcf.gz \
        -L ${interval}
    """
}


// ----------------------------
// GENOMICSDB IMPORT
// ----------------------------

process WES_GENOMICSDB_IMPORT {

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
    gatk GenomicsDBImport \
        --genomicsdb-workspace-path ${interval_name}_db \
        -R ${fasta} \
        -L ${intervals} \
        ${gvcfs.collect { "-V ${it}" }.join(' ')}
    """
}


process WGS_GENOMICSDB_IMPORT {

    tag "${interval_name}"
    label 'gatk'

    input:
        tuple val(interval_name), path(gvcfs)
        tuple val(interval_name), path(tbis)
        tuple path(fasta), path(fai), path(fasta_dict)

    output:
        tuple val(interval_name), path("${interval_name}_db")

    script:
    """
    gatk GenomicsDBImport \
        --genomicsdb-workspace-path ${interval_name}_db \
        -R ${fasta} \
        ${gvcfs.collect { "-V ${it}" }.join(' ')}
    """
}


// ----------------------------
// GENOTYPING
// ----------------------------

process WES_GENOTYPE_GVCF {

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
    gatk GenotypeGVCFs \
        -R ${fasta} \
        -V gendb://${gendb} \
        -O ${interval_name}.vcf.gz \
        -L ${intervals} \
        --merge-input-intervals
    """
}


process WGS_GENOTYPE_GVCF {

    tag "${interval_name}"
    label 'gatk'

    input:
        tuple val(interval_name), path(gendb)
        tuple path(fasta), path(fai), path(fasta_dict)

    output:
        path("${interval_name}.vcf.gz")

    script:
    """
    gatk GenotypeGVCFs \
        -R ${fasta} \
        -V gendb://${gendb} \
        -O ${interval_name}.vcf.gz \
        -L ${interval_name} \
        --merge-input-intervals
    """
}


// ----------------------------
// POSTPROCESSING
// ----------------------------

process GATHER_VCF {

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


process VARIANT_FILTERING {

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
    bcftools norm -a -m - -f ${fasta} -Oz -o cohort_hc_norm.vcf.gz

    tabix -p vcf cohort_hc_norm.vcf.gz
    """
}


process HAPLOTYPE_CALLER_EXTRACT_GT {

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