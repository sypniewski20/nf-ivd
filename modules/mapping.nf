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

