process SPLICE_AI {
    publishDir "${params.outfolder}/${params.runID}/spliceAI/", mode: 'copy', overwrite: true
    label 'spliceAI'
	label 'mem_36GB'
	label 'core_36'
	input:
		path(vcf)
        path(tbi)
        tuple path(fasta), path(fai)

	output:
		tuple path("${vcf.baseName}_spliceAI.vcf.gz"), path("${vcf.baseName}_spliceAI.vcf.gz.tbi")
	script:
      def genome = (fasta.name =~ /(?i)GRCh38|hg38|Homo_sapiens_assembly38/) ? "grch38" : "grch37"
      def distance = params.spliceai_distance ?: 50
      
		"""

        spliceai -I ${vcf} \
                 -O ${vcf.baseName}_spliceAI.vcf \
                 -R ${fasta} \
                 -A ${genome} \
                 -D ${distance}


        bgzip ${vcf.baseName}_spliceAI.vcf
        tabix -p vcf ${vcf.baseName}_spliceAI.vcf.gz

		"""

}

process VEP_GERMLINE_SNV {
    publishDir "${params.outfolder}/${params.runID}/vep/", mode: 'copy', overwrite: true
    label 'vep'
    label 'mem_36GB'
    label 'core_36'

    input:
        path(vcf)
        path(tbi)
        tuple path(spliceai_vcf), path(spliceai_tbi)
        tuple path(dbnsfp), path(dbnsfp_tbi)
        tuple path(fasta), path(fai)

    output:
        path("${vcf.baseName}.vep.tsv.gz")

    script:
        def genome  = (fasta.name =~ /(?i)GRCh38|hg38|Homo_sapiens_assembly38/) ? "GRCh38" : "GRCh37"
        def clinvar = params.clinvar
            ? "--custom file=${params.clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%MC%CLNDISDB%CLNDISDBINC"
            : ""
        def dbNSFP  = params.dbnsfp
            ? "--plugin dbNSFP,${params.dbnsfp},,SIFT_score,SIFT_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,MetaRNN_score,MetaRNN_pred,REVEL_score,BayesDel_addAF_score,BayesDel_addAF_pred,BayesDel_noAF_score,BayesDel_noAF_pred,ClinPred_score,ClinPred_pred,AlphaMissense_score,AlphaMissense_pred,CADD_raw,CADD_phred,GERP++_RS,phyloP100way_vertebrate,gnomAD4.1_joint_NFE_AF,gnomAD4.1_joint_AF"
            : ""
        def spliceai = (params.seq_type != 'WES' && spliceai_vcf.name != 'NO_FILE_VCF')
            ? "--custom file=${spliceai_vcf},short_name=SpliceAI,format=vcf,type=overlap,coords=0,fields=ALLELE%SYMBOL%DS_AG%DS_AL%DS_DG%DS_DL"
            : ""
        """
        vep \
            --cache \
            --dir_cache ${params.vep_cache} \
            --species homo_sapiens \
            --assembly ${genome} \
            --buffer_size 50000 \
            --input_file ${vcf} \
            --output_file "${vcf.baseName}.vep.tsv.gz" \
            --format vcf \
            --compress_output bgzip \
            --fasta ${fasta} \
            --tab \
            --no_stats \
            --fork ${task.cpus} \
            --force_overwrite \
            --hgvs \
            --hgvsg \
            --symbol \
            --numbers \
            --domains \
            --protein \
            --biotype \
            --uniprot \
            --variant_class \
            ${spliceai} \
            ${clinvar} \
            ${dbNSFP} \
            --offline \
            --cache_version ${params.cache_version} \
            --pick \
            --pick_order mane_select,mane_plus_clinical,canonical,tsl,biotype,rank
        """
}

process FILTER_VEP {
    publishDir "${params.outfolder}/${params.runID}/vep/", mode: 'copy', overwrite: true
    label 'vep'
	label 'mem_36GB'
	label 'core_36'
	input:
		path(tsv)
	output:
		path("${tsv.baseName}.vep.filtered.tsv.gz")
	script:
		"""

        filter_vep \
            -i ${tsv} \
            --format tab \
            --filter "PICK == 1" \
            -o - | bgzip -c > ${tsv.baseName}.vep.filtered.tsv.gz
        
		"""

}