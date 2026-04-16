process DRAGMAP_BAM {
	publishDir "${params.outfolder}/${params.runID}/BAM/", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'large'
	input:
		tuple val(sample), val(LB), path(read_1), path(read_2)
		path(fasta_dir)
		path(fasta)
	output:
		tuple val(sample), path("${sample}_sorted.bam"), path("${sample}_sorted.bam.bai"), emit: ch_bam
		tuple val(sample), path("${sample}_sorted.flagstat"), emit: ch_flagstat
		tuple val(sample), path("${sample}_sorted.bam.md5"), emit: ch_md5
	script:
		"""

		dragen-os \
			--num-threads ${task.cpus} \
			-r ${fasta_dir} \
			-1 ${read_1} \
			-2 ${read_2} \
			--RGSM ${sample} \
			--RGID ${sample} | \
		samtools view --reference ${fasta} \
					  --threads ${task.cpus} -b | \
		samtools sort -@ ${task.cpus} \
					  -O bam > ${sample}_sorted.bam

		samtools index --nthreads ${task.cpus} ${sample}_sorted.bam

		samtools quickcheck ${sample}_sorted.bam
		md5sum ${sample}_sorted.bam > ${sample}_sorted.bam.md5
		samtools flagstat ${sample}_sorted.bam > ${sample}_sorted.flagstat

		"""
}