process EXCLUSION_BED {
    label "core"
    input:
        path bed
        tuple path(fasta), path(fai)
    output:
        path "exclusion.bed"
    script:
    """
    
    bedtools sort -i ${bed} -g ${fai} \
        | bedtools complement -i - -g ${fai} \
        > exclusion.bed
    """
}

process DELLY_CALL {
    tag "${sample}"
    label 'delly'
    input:
        tuple val(sample), path(bam), path(bai)
        tuple path(fasta), path(fai)
        path exclusion_bed
    output:
        tuple val(sample), path("${sample}_delly.bcf"), path("${sample}_delly.bcf.csi"), emit: bcf
    script:
    """
    delly call \
        -g ${fasta} \
        -x ${exclusion_bed} \
        -o ${sample}_delly.bcf \
        ${bam}

    """
}

process DELLY_SITE_MERGE {
    label 'delly'
    input:
        path bcfs
        path csis
    output:
        path "delly_merged_sites.bcf", emit: sites
    script:
    """
    delly merge \
        -o delly_merged_sites.bcf \
        ${bcfs}
    """
}

process DELLY_GENOTYPE {
    tag "${sample}"
    label 'delly'
    input:
        tuple val(sample), path(bam), path(bai)
        tuple path(fasta), path(fai)
        path exclusion_bed
        path merged_sites
    output:
        tuple val(sample), path("${sample}_delly_genotyped.bcf"), path("${sample}_delly_genotyped.bcf.csi"), emit: bcf
    script:
    """
    delly call \
        -g ${fasta} \
        -x ${exclusion_bed} \
        -v ${merged_sites} \
        -o ${sample}_delly_genotyped.bcf \
        ${bam}

    """
}

process DELLY_FILTER {
    label 'delly'
    input:
        path bcfs  // all per-sample genotyped BCFs collected
        path csis
    output:
        path "delly_filtered.bcf",      emit: bcf
        path "delly_filtered.bcf.csi",  emit: csi
    script:
    """
    bcftools merge \
        --merge none \
        --force-single \
        -Ob \
        -o delly_cohort.bcf \
        ${bcfs}
    bcftools index delly_cohort.bcf

    delly filter \
        -f germline \
        -o delly_filtered.bcf \
        delly_cohort.bcf

    """
}

process DELLY_TO_VCF {
    publishDir "${params.outfolder}/${params.runID}/SV/delly", mode: 'copy', overwrite: true
    label 'delly'
    input:
        path bcf
        path csi
    output:
        tuple path("delly_cohort.vcf.gz"), path("delly_cohort.vcf.gz.tbi"), emit: vcf
    script:
    """
    bcftools view \
        -Oz \
        -o delly_cohort.vcf.gz \
        ${bcf}

    tabix -p vcf delly_cohort.vcf.gz
    """
}

process DELLY_TO_TSV {
    publishDir "${params.outfolder}/${params.runID}/SV/delly", mode: 'copy', overwrite: true
    label 'delly'
    input:
        tuple path(vcf), path(tbi)
    output:
        path "delly_cohort.tsv"
    script:
    """
    
    echo -e "CHROM\tSTART\tEND\tID\tSVTYPE\tSVLEN\tFILTER\tCIPOS\tCIEND\tSAMPLE\tGT\tDR\tDV\tRR\tRV" > delly_cohort.tsv
    bcftools query \
        -f '[%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%SVLEN\t%FILTER\t%CIPOS\t%CIEND\t%SAMPLE\t%GT\t%DR\t%DV\t%RR\t%RV\n]' \
        ${vcf} >> delly_cohort.tsv
    

    """
}