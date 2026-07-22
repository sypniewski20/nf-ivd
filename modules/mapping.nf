process DRAGMAP_BAM {
    tag "${sample}"
    label 'core'
    label 'xlarge'
    input:
        // tuple contains: sample name, Library ID (LB), Platform (PL), and FASTQ paths
        tuple val(sample), val(ID), val(LB), val(PL), val(PU), path(read_1), path(read_2)
        tuple path(fasta_dir), path(fasta), path(fasta_fai)

    output:
        tuple val(sample), path("${sample}_sorted.bam"), path("${sample}_sorted.bam.bai")

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
        
        """
}

process BWA_BAM {
    tag "${sample}"
    label 'core'
    label 'xlarge'
    input:
        tuple val(sample), val(ID), val(LB), val(PL), val(PU), path(read_1), path(read_2)
        tuple path(fasta), path(fasta_fai), path(fasta_amb), path(fasta_ann), path(fasta_bwt), path(fasta_sa), path(fasta_pac)

    output:
        tuple val(sample), path("${sample}_sorted.bam"), path("${sample}_sorted.bam.bai")

    script:
        def rg_id = (ID && ID.trim()) ? ID : [sample, LB, PU].findAll { it && it.trim() }.join('_')
        def rg_line = "@RG\\tID:${rg_id}\\tSM:${sample}\\tLB:${LB}\\tPL:${PL}\\tPU:${PU}"
        """
        #!/bin/bash
        set -eo pipefail
        
		bwa mem -t ${task.cpus} \
			-R "${rg_line}" \
			${fasta} \
			${read_1} \
			${read_2} | \
		samtools view --reference ${fasta} \
					  --threads ${task.cpus} -b | \
		samtools sort -@ ${task.cpus} \
					  -O bam > ${sample}_sorted.bam

        samtools index -@ ${task.cpus} ${sample}_sorted.bam
        
        """
}

process MARK_DUPLICATES {
    publishDir "${params.outfolder}/${params.runID}/BAM",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            "${sample}/${filename}"
        }    
    tag "${sample}"
    label 'core'
    label 'xlarge'
    input:
        tuple val(sample), path(bam), path(bai)

    output:
        tuple val(sample), path("${sample}_sorted_markdup.bam"), path("${sample}_sorted_markdup.bam.bai"), emit: ch_bam
        tuple path("${sample}_sorted_markdup.flagstat"), path("${sample}_sorted_markdup.metrics"), emit: ch_metrics
        tuple val(sample), path("${sample}_sorted_markdup.bam.md5"), emit: ch_md5

    script:
        """
        #!/bin/bash
        set -eo pipefail
        
        gatk MarkDuplicates \
            -I ${bam} \
            -O ${sample}_sorted_markdup.bam \
            -M ${sample}_sorted_markdup.metrics \
            --TMP_DIR . 

        # Clinical Integrity Checks

        samtools index ${sample}_sorted_markdup.bam

        samtools quickcheck ${sample}_sorted_markdup.bam
        md5sum ${sample}_sorted_markdup.bam > ${sample}_sorted_markdup.bam.md5
        samtools flagstat -@ ${task.cpus} ${sample}_sorted_markdup.bam > ${sample}_sorted_markdup.flagstat
        
        """
}

