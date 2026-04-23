process DRAGMAP_BAM {
    publishDir "${params.outfolder}/${params.runID}/BAM/", mode: 'copy', overwrite: true
    tag "${sample}"
    label 'gatk'
    label 'xlarge'
    input:
        // tuple contains: sample name, Library ID (LB), Platform (PL), and FASTQ paths
        tuple val(sample), val(ID), val(LB), val(PL), val(PU), path(read_1), path(read_2)
        tuple path(fasta_dir), path(fasta), path(fasta_fai)

    output:
        tuple val(sample), path("${sample}_sorted.bam"), path("${sample}_sorted.bam.bai"), emit: ch_bam
        tuple val(sample), path("${sample}_sorted.flagstat"),                            emit: ch_flagstat
        tuple val(sample), path("${sample}_sorted.bam.md5"),                            emit: ch_md5

    script:
        def rg_id = (ID && ID.trim()) ? ID : [sample, LB, PU].findAll { it && it.trim() }.join('_')
        def rg_line = "@RG\\tID:${rg_id}\\tSM:${sample}\\tLB:${LB}\\tPL:${PL}\\tPU:${PU}"
        """
        #!/bin/bash
        set -eo pipefail
        
		dragen-os \
			-r ${fasta_dir} \
			-1 ${read_1} \
			-2 ${read_2} \
            --num-threads ${task.cpus} | \
        samtools addreplacerg \
            -r "${rg_line}" \
            -m overwrite_all \
            -O bam - | \
		samtools sort -@ ${task.cpus} \
					  -O bam \
                      -o ${sample}_sorted.bam

		samtools index -@ ${task.cpus} ${sample}_sorted.bam

        # Clinical Integrity Checks
        samtools quickcheck ${sample}_sorted.bam
        md5sum ${sample}_sorted.bam > ${sample}_sorted.bam.md5
        samtools flagstat -@ ${task.cpus} ${sample}_sorted.bam > ${sample}_sorted.flagstat
        
        """
}

process DRAGMAP_STREAM_BAM {
    tag "${sample}_${LB}"
    label 'gatk'
    label 'xlarge'
    input:
        // tuple contains: sample name, Library ID (LB), Platform (PL), and FASTQ paths
        tuple val(sample), val(ID), val(LB), val(PL), val(PU), path(R1_URL), path(R2_URL)
        tuple path(fasta_dir), path(fasta), path(fasta_fai)

    output:
        tuple val(sample), path("${sample}_${LB}_sorted.bam"), path("${sample}_${LB}_sorted.bam.bai")

    script:
        def rg_id = (ID && ID.trim()) ? ID : [sample, LB, PU].findAll { it && it.trim() }.join('_')
        def rg_line = "@RG\\tID:${rg_id}\\tSM:${sample}\\tLB:${LB}\\tPL:${PL}\\tPU:${PU}"
        """
        #!/bin/bash
        set -eo pipefail
        
		dragen-os \
			-r ${fasta_dir} \
			-1 <(curl -sL "${R1_URL}") \
			-2 <(curl -sL "${R2_URL}") \
            --num-threads ${task.cpus} | \
        samtools addreplacerg \
            -r "${rg_line}" \
            -m overwrite_all \
            -O bam - | \
		samtools sort -@ ${task.cpus} \
					  -O bam \
                      -o ${sample}_${LB}_sorted.bam

		samtools index -@ ${task.cpus} ${sample}_${LB}_sorted.bam

        # Clinical Integrity Checks
        samtools quickcheck ${sample}_${LB}_sorted.bam

        
        """
}

process MERGE_BAM {
    publishDir "${params.outfolder}/${params.runID}/BAM/", mode: 'copy', overwrite: true
    tag "${sample}"
    label 'gatk'
    label 'xlarge'
    input:
        tuple val(sample), val(lb), path(bams)
        tuple val(sample), val(lb), path(bais)

    output:
        tuple val(sample), path("${sample}_merged.bam"), path("${sample}_merged.bam.bai"), emit: ch_bam
        tuple val(sample), path("${sample}_merged.flagstat"),                              emit: ch_flagstat
        tuple val(sample), path("${sample}_merged.bam.md5"),                              emit: ch_md5

    script:
        """
        #!/bin/bash
        set -eo pipefail

        samtools merge \
            --threads ${task.cpus} \
            -f ${sample}_merged.bam \
            ${bams}

        samtools index -@ ${task.cpus} ${sample}_merged.bam

        samtools quickcheck ${sample}_merged.bam
        md5sum ${sample}_merged.bam > ${sample}_merged.bam.md5
        samtools flagstat -@ ${task.cpus} ${sample}_merged.bam > ${sample}_merged.flagstat
        """
}