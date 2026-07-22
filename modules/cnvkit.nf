// modules/cnvkit.nf

process CNVKIT_TARGET {
    label 'cnvkit'

    input:
        path(bed)

    output:
        path("targets.bed")

    script:
        """
        cnvkit.py target ${bed} --split -o targets.bed
        """
}

process CNVKIT_ANTITARGET {
    label 'cnvkit'

    input:
        path(target_bed)

    output:
        path("antitargets.bed")

    script:
        """
        cnvkit.py antitarget ${target_bed} -o antitargets.bed
        """
}

process CNVKIT_COVERAGE {
    tag "${sample}"
    label 'cnvkit'
    label 'mem_8GB'
    label 'core_8'

    input:
        tuple val(sample), path(bam), path(bai)
        path(targets)
        path(antitargets)

    output:
        tuple val(sample), path("${sample}.targetcoverage.cnn"), path("${sample}.antitargetcoverage.cnn")

    script:
        """
        cnvkit.py coverage ${bam} ${targets} \
            -o ${sample}.targetcoverage.cnn \
            -p ${task.cpus}

        cnvkit.py coverage ${bam} ${antitargets} \
            -o ${sample}.antitargetcoverage.cnn \
            -p ${task.cpus}
        """
}

process CNVKIT_REFERENCE {
    label 'cnvkit'
    label 'mem_8GB'

    input:
        tuple path(fasta), path(fai)
        path(targets)
        path(antitargets)

    output:
        path("pooled_reference.cnn")

    script:
        """
        cnvkit.py reference \
            -t ${targets} \
            -a ${antitargets} \
            -f ${fasta} \
            -o pooled_reference.cnn
        """
}

process CNVKIT_FIX {
    tag "${sample}"
    label 'cnvkit'

    input:
        tuple val(sample), path(target_cnn), path(antitarget_cnn)
        path(reference)

    output:
        tuple val(sample), path("${sample}.cnr"), emit: cnr

    script:
        def edge_flag = (params.seq_type == 'WES') ? "" : "--no-edge"
        """
        cnvkit.py fix \
            ${target_cnn} \
            ${antitarget_cnn} \
            ${reference} \
            ${edge_flag} \
            -o ${sample}.cnr
        """
}

process CNVKIT_SEGMENT {
    tag "${sample}"
    label 'cnvkit'
    label 'mem_8GB'
    label 'core_4'

    input:
        tuple val(sample), path(cnr)

    output:
        tuple val(sample), path("${sample}.cns"), emit: cns

    script:
        """
        cnvkit.py segment ${cnr} \
            -o ${sample}.cns \
            -p ${task.cpus} \
            --smooth-cbs
        """
}

process CNVKIT_CALL {
    tag "${sample}"
    label 'cnvkit'
    publishDir path: { "${params.outfolder}/${params.runID}/SV/cnvkit/${sample}" }, mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(cns)

    output:
        tuple val(sample), path("${sample}.call.cns"), emit: cns

    script:
        """
        cnvkit.py call ${cns} \
            --ploidy 2 \
            --method clonal \
            -o ${sample}.call.cns
        """
}

process CNVKIT_SCATTER {
    tag "${sample}"
    label 'cnvkit'
    publishDir path: { "${params.outfolder}/${params.runID}/SV/cnvkit/${sample}" }, mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(cnr), path(cns)

    output:
        tuple val(sample), path("${sample}_scatter.pdf"), emit: pdf

    script:
        """
        cnvkit.py scatter ${cnr} \
            -s ${cns} \
            -o ${sample}_scatter.pdf
        """
}

process CNVKIT_DIAGRAM {
    tag "${sample}"
    label 'cnvkit'
    publishDir path: { "${params.outfolder}/${params.runID}/SV/cnvkit/${sample}" }, mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(cnr), path(cns)

    output:
        tuple val(sample), path("${sample}_diagram.pdf"), emit: pdf

    script:
        """
        cnvkit.py diagram ${cnr} \
            -s ${cns} \
            -o ${sample}_diagram.pdf
        """
}

process CNVKIT_TO_VCF {
    tag "${sample}"
    label 'cnvkit'
    publishDir path: { "${params.outfolder}/${params.runID}/SV/cnvkit/${sample}" }, mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(cns)

    output:
        tuple path("${sample}_cnvkit.vcf.gz"), path("${sample}_cnvkit.vcf.gz.tbi")

    script:
        """
        cnvkit.py export vcf ${cns} -o ${sample}_cnvkit.vcf

        bgzip ${sample}_cnvkit.vcf
        tabix -p vcf ${sample}_cnvkit.vcf.gz
        """
}

process MERGE_VCF {
    label 'cnvkit'
    publishDir "${params.outfolder}/${params.runID}/SV/cnvkit", mode: 'copy', overwrite: true

    input:
        path(vcf)
        path(tbi)

    output:
        tuple path("cnvkit_merged.vcf.gz"), path("cnvkit_merged.vcf.gz.tbi")

    script:
        """
        bcftools merge ${vcf} -Ou | \
        bcftools annotate --set-id +'%CHROM\\_%POS\\_%SVTYPE' -Ou | \
        bcftools +fill-tags -Ou -- -t AF,AC | \
        bcftools sort -Oz -o cnvkit_merged.vcf.gz

        tabix -p vcf cnvkit_merged.vcf.gz
        """
}