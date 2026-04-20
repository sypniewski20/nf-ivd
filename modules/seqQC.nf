process FASTP_PROCESSING {
	publishDir "${params.outfolder}/${params.runID}/fastp/${sample}", pattern: "fastp.*", mode: 'copy', overwrite: true
	label 'gatk'
	tag "${sample}"
	label 'medium'
	input:
		tuple val(sample), val(ID), val(LB), val(PL), val(PU), path(read_1), path(read_2)
	output:
		tuple val(sample), val(ID), val(LB), val(PL), val(PU), path("${sample}.filtered.R1.fq.gz"), path("${sample}.filtered.R2.fq.gz"), emit: fastq_filtered
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

process FASTP_STREAM {
	publishDir "${params.outfolder}/${params.runID}/fastp/${sample}", pattern: "fastp.*", mode: 'copy', overwrite: true
    tag "${sample}_${ID}"
	label 'gatk'
    label 'tiny'
    input:
		tuple val(sample), val(ID), val(LB), val(PL), val(PU), path(R1_URL), path(R2_URL)
    output:
        tuple val(sample), val(ID), val(LB), val(PL), val(PU), path("${sample}_${ID}_R1_fastp.fq.gz"), path("${sample}_${ID}_R2_fastp.fq.gz"), emit: fastq_filtered
		tuple path("${sample}_${ID}_fastp.html"), path("${sample}_${ID}_fastp.json"), emit: fastp_log
    script:
        """

        fastp -i <(curl -sL "${R1_URL}") \
              -I <(curl -sL "${R2_URL}") \
              -o ${sample}_${ID}_R1_fastp.fq.gz \
              -O ${sample}_${ID}_R2_fastp.fq.gz \
              -w ${task.cpus} \
			  --html ${sample}_${ID}_fastp.html \
			  --json ${sample}_${ID}_fastp.json \
              --detect_adapter_for_pe
			  
        """
}

process MOSDEPTH {
	publishDir "${params.outfolder}/${params.runID}/BAMQC", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'qc'
	label 'medium'
	input:
		tuple val(sample), path(bam), path(bai)
		tuple path(fasta), path(fai)
	output:
		path("*")
	script:
		"""

		mosdepth -f ${fasta} -n --fast-mode --by 500  ${sample} ${bam} --threshold 1,10,20,30 -t ${task.cpus}

		"""
}

process MOSDEPTH_EXOME {
	publishDir "${params.outfolder}/${params.runID}/BAMQC", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'qc'
	label 'medium'
	input:
		tuple val(sample), path(bam), path(bai)
		path(intervals)
		tuple path(fasta), path(fai)
	output:
		path("*")
	script:
		"""

		mosdepth -f ${fasta} -n --threshold 1,10,20,30 --by ${intervals} --fast-mode ${sample} ${bam} -t ${task.cpus}

		"""
}

process FASTQC {
	publishDir "${params.outfolder}/${params.runID}/fastqc/${sample}", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'qc'
	label 'small'
	input:
		tuple val(sample), val(ID), val(LB), val(PL), val(PU), path(read_1), path(read_2)
	output:
		path("*")
	script:
		"""

		fastqc ${read_1} ${read_2} -t ${task.cpus}

		"""
}

process MULTIQC {
	publishDir "${params.outfolder}/${params.runID}/multiqc", mode: 'copy', overwrite: true
	label 'qc'
    label 'tiny'
	input:
		path(fastqc)
		path(flagstat)
		path(mosdepth)
		path(vcf_stats)
	output:
		path("${params.runID}*")
	script:
		"""

		multiqc . --filename ${params.runID} --verbose

		"""
}