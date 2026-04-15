// ============================================================
// Tumor-only Mutect2 module (WITH PON, WES/WGS VERSION)
// ============================================================


// ----------------------------
// MUTECT2 CALLING
// ----------------------------

process MUTECT2_SOMATIC_ONLY {

    label 'gatk'
    label 'mem_4GB'
    label 'core_8'
    tag "${sample}"

    input:
        tuple val(sample), path(bam), path(bai)
        tuple path(fasta), path(fai), path(fasta_dict)
        tuple path(snv_resource), path(snv_tbi)
        tuple path(pon), path(pon_tbi)              // ✅ ADDED PON
        val(interval_padding)
        path(intervals)

    output:
        tuple val(sample),
              path("${sample}.mutect2.vcf.gz"),
              path("${sample}.mutect2.vcf.gz.tbi"),
              path("${sample}.f1r2.tar.gz"),
              emit: vcf

    script:
    """
    gatk Mutect2 \
        -R ${fasta} \
        -I ${bam} \
        -L ${intervals} \
        --interval-padding ${interval_padding} \
        --germline-resource ${snv_resource} \
        --panel-of-normals ${pon} \          # ✅ ADDED
        -O ${sample}.mutect2.vcf.gz \
        --native-pair-hmm-threads ${task.cpus} \
        --max-mnp-distance 0 \
        --dont-use-soft-clipped-bases true \
        --f1r2-tar-gz ${sample}.f1r2.tar.gz
    """
}


// ----------------------------
// ORIENTATION MODEL
// ----------------------------

process MUTECT2_ORIENTATION_MODEL {

    label 'gatk'
    label 'mem_4GB'
    label 'core_8'
    tag "${sample}"

    input:
        tuple val(sample), path(f1r2)

    output:
        tuple val(sample), path("${sample}.read-orientation-model.tar.gz")

    script:
    """
    gatk LearnReadOrientationModel \
        -I ${f1r2} \
        -O ${sample}.read-orientation-model.tar.gz
    """
}


// ----------------------------
// ARTIFACT METRICS
// ----------------------------

process MUTECT2_ARTIFACT_METRICS {

    publishDir "${params.outfolder}/${params.runID}/SNV/mutect2/artifacts",
        mode: 'copy', overwrite: true

    label 'gatk'
    label 'mem_8GB'
    label 'core_8'
    tag "${sample}"

    input:
        tuple val(sample), path(bam), path(bai)
        tuple path(fasta), path(fai), path(fasta_dict)

    output:
        path("${sample}.artifact-metrics*")

    script:
    """
    gatk CollectSequencingArtifactMetrics \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample}.artifact-metrics
    """
}


// ----------------------------
// FILTERING
// ----------------------------

process MUTECT2_FILTER {

    publishDir "${params.outfolder}/${params.runID}/SNV/mutect2/",
        mode: 'copy', overwrite: true

    label 'gatk'
    label 'mem_4GB'
    label 'core_8'
    tag "${sample}"

    input:
        tuple val(sample), path(vcf), path(tbi), path(rom)
        tuple path(fasta), path(fai), path(fasta_dict)

    output:
        tuple val(sample),
              path("${sample}.PASS.vcf.gz"),
              path("${sample}.PASS.vcf.gz.tbi")

    script:
    """
    gatk FilterMutectCalls \
        -V ${vcf} \
        -R ${fasta} \
        --ob-priors ${rom} \
        -O ${sample}.filtered.vcf.gz

    bcftools annotate \
        --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' \
        ${sample}.filtered.vcf.gz \
        -Oz -o ${sample}.PASS.vcf.gz

    tabix -fp vcf ${sample}.PASS.vcf.gz
    """
}


// ----------------------------
// GENOTYPE EXTRACTION
// ----------------------------

process MUTECT2_EXTRACT_GT {

    tag "${sample}"

    publishDir "${params.outfolder}/${params.runID}/SNV/mutect2/",
        mode: 'copy', overwrite: true

    label 'gatk'
    label 'mem_16GB'
    label 'core_4'

    input:
        tuple val(sample), path(vcf), path(tbi)

    output:
        path("${sample}.gt.tsv.gz")

    script:
    """
    echo -e "ID\\tSAMPLE\\tFILTER\\tDP\\tAF\\tGT" | \
    bgzip -c > ${sample}.gt.tsv.gz

    bcftools query -f "[%ID\\t%SAMPLE\\t%FILTER\\t%DP\\t%AF\\t%GT\\n]" ${vcf} | \
    bgzip -c >> ${sample}.gt.tsv.gz
    """
}