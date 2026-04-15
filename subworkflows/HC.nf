include {
    WES_GVCF_HAPLOTYPE_CALLER;  WGS_GVCF_HAPLOTYPE_CALLER
    WES_HAPLOTYPE_CALLER;       WGS_HAPLOTYPE_CALLER
    WES_GENOTYPE_GVCF;          WGS_GENOTYPE_GVCF
    WES_GENOMICSDB_IMPORT;      WGS_GENOMICSDB_IMPORT
    HAPLOTYPE_CALLER_EXTRACT_GT
    GATHER_VCF
    VARIANT_FILTERING
} from "../modules/HaplotypeCaller.nf"
workflow germline_workflow {

    take:
        ch_bam
        params

    main:

        def isWES    = params.seq_type == 'WES'
        def isCohort = params.calling_mode == 'cohort'

        def gvcf = isWES
            ? WES_GVCF_HAPLOTYPE_CALLER(ch_bam, params.fasta, params.interval_padding)
            : WGS_GVCF_HAPLOTYPE_CALLER(ch_bam, params.fasta)

        def db = isWES
            ? WES_GENOMICSDB_IMPORT(gvcf.vcf, gvcf.out.tbi, params.fasta, params.bed)
            : WGS_GENOMICSDB_IMPORT(gvcf.vcf, gvcf.out.tbi, params.fasta)

        def genotype = isWES
            ? WES_GENOTYPE_GVCF(db, params.fasta, params.bed)
            : WGS_GENOTYPE_GVCF(db, params.fasta)

        def hc = isWES
            ? WES_HAPLOTYPE_CALLER(ch_bam, params.fasta, params.interval_padding)
            : WGS_HAPLOTYPE_CALLER(ch_bam, params.fasta)

        GATHER_VCF(genotype.mix(hc.vcf).collect())

        VARIANT_FILTERING(GATHER_VCF.out, params.fasta)

        HAPLOTYPE_CALLER_EXTRACT_GT(VARIANT_FILTERING.out)

    emit:
        hc_vcf = VARIANT_FILTERING.out
}