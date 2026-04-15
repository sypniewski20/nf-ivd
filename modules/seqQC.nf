process FASTP_PROCESSING {
	publishDir "${params.outfolder}/${params.runID}/fastp/${sample}", pattern: "fastp.*", mode: 'copy', overwrite: true
	label 'qc'
	tag "${sample}"
	label 'mem_8GB'
	label 'core_8'
	input:
		tuple val(sample), val(LB), path(read_1), path(read_2)
	output:
		tuple val(sample), val(LB), path("${sample}.filtered.R1.fq.gz"), path("${sample}.filtered.R2.fq.gz"), emit: fastq_filtered
		tuple path("${sample}_fastp.html"), path("${sample}_fastp.json"), emit: fastp_log
	script:
		"""
		fastp -i ${read_1} \
			  -I ${read_2} \
			  -o ${sample}.filtered.R1.fq.gz \
			  -O ${sample}.filtered.R2.fq.gz \
			  -t ${task.cpus} \
			  --html ${sample}_fastp.html \
			  --json ${sample}_fastp.json \
	  		  --detect_adapter_for_pe
		"""
}

process MOSDEPTH {
	publishDir "${params.outfolder}/${params.runID}/BAMQC", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'qc'
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
	output:
		path("*")
	script:
		"""

		mosdepth -f ${fasta} -n --fast-mode --by 500  ${sample} ${bam} --threshold 1,10,20,30 -t ${task.cpus}
		python /mosdepth/plot-dist.py ${sample}.mosdepth.global.dist.txt
		mv dist.html ${sample}_mosdepth.html

		"""
}

process MOSDEPTH_EXOME {
	publishDir "${params.outfolder}/${params.runID}/BAMQC", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'qc'
	label 'mem_16GB'
	label 'core_4'
	input:
		tuple val(sample), path(bam), path(bai)
		path(intervals)
		path(fasta)
	output:
		path("*")
	script:
		"""

		mosdepth -f ${fasta} -n --threshold 1,10,20,30 --by ${intervals} --fast-mode ${sample} ${bam} -t ${task.cpus}
		python /mosdepth/plot-dist.py ${sample}.mosdepth.region.dist.txt
		mv dist.html ${sample}_mosdepth.html

		"""
}

process FASTQC {
	publishDir "${params.outfolder}/${params.runID}/fastqc/${sample}", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_8GB'
	label 'core_1'
	input:
		tuple val(sample), val(LB), path(read_1), path(read_2)
	output:
		path("*")
	script:
		"""

		fastqc ${read_1} ${read_2} -t ${task.cpus}

		"""
}

process MULTIQC {
	publishDir "${params.outfolder}/${params.runID}/multiqc", mode: 'copy', overwrite: true
	label 'gatk'
	label 'mem_8GB'
	label 'core_1'
	input:
		path(reports)
	output:
		path("${params.runID}*")
	script:
		"""

		multiqc . --filename ${params.runID} --verbose

		"""
}