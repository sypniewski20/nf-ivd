process DEEP_VARIANT {
    tag "${sample}"
    label params.gpus == 0 ? 'deepvariant' : 'deepvariant_gpu'
    label 'xlarge'
    publishDir "${params.outfolder}/${params.runID}/deepvariant", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(bam), path(bai)
        tuple path(fasta), path(fai)

    output:
        tuple val(sample), path("${sample}_dv.vcf.gz"), path("${sample}_dv.vcf.gz.tbi"), emit: vcf
        tuple val(sample), path("${sample}_dv.gvcf.gz"), path("${sample}_dv.gvcf.gz.tbi"), emit: gvcf

    script:
        """
        /opt/deepvariant/bin/run_deepvariant \
            --model_type ${params.seq_type} \
            --ref ${fasta} \
            --reads ${bam} \
            --output_vcf ${sample}_dv.vcf.gz \
            --output_gvcf ${sample}_dv.gvcf.gz \
            --num_shards ${task.cpus}
        """
}

process GLNEXUS {
    label 'core'
    label 'large'
    publishDir "${params.outfolder}/${params.runID}/deepvariant", mode: 'copy', overwrite: true

    input:
        path(gvcf)
        path(gvcf_tbi)
        tuple path(fasta), path(fai)

    output:
        path("norm_glnexus_deepvariant.vcf.gz"), emit: vcf
        path("norm_glnexus_deepvariant.vcf.gz.tbi"), emit: tbi

    script:
        """
        LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so

        glnexus_cli \
        --threads ${task.cpus} \
        --config DeepVariant${params.seq_type} \
        ${gvcf} | \
        bcftools norm -a --atom-overlaps . -m - -f ${fasta} -Ou | \
        bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' -Ou | \
        bcftools +fill-tags -Ou -- -t AF,AC | \
        bcftools sort -Oz -o norm_glnexus_deepvariant.vcf.gz

        tabix -p vcf norm_glnexus_deepvariant.vcf.gz
        """
}

process DV_EXTRACT_GT {
    label 'tiny'
    label 'core'
    publishDir "${params.outfolder}/${params.runID}/deepvariant", mode: 'copy', overwrite: true

    input:
        path(vcf)
        path(tbi)

    output:
        path("${vcf.simpleName}_gt.table"), emit: ch_gt_table
        path("${vcf.simpleName}_gt.table.md5"), emit: ch_gt_table_md5

    script:
        """
        gatk VariantsToTable \
            -V ${vcf} \
            -F CHROM -F POS -F ID -F REF -F ALT -GF GT -GF DP \
            -O ${vcf.simpleName}_gt.table

        md5sum ${vcf.simpleName}_gt.table > ${vcf.simpleName}_gt.table.md5
        """
}