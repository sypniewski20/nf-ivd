process MERGE_SVS {
    label 'gatk'
    label 'small'
    input:
        path(manta_vcf)
        path(gcnv_vcf)

    output:
        path("merged_consensus.vcf")

    script:
    """

    echo "${manta_vcf}" > sample_list.txt
    echo "${gcnv_vcf}" >> sample_list.txt

    SURVIVOR merge sample_list.txt 500 2 1 0 50 merged_consensus.vcf
    """
}