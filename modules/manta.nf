process JOINT_DIPLOID_ANALYSIS {
    label 'medium'
    label 'manta'
    publishDir "${params.outfolder}/${params.runID}/", mode: 'copy'
    input:
        path bam_files
        path bai_files
        path fasta

    output:
        path "manta"

    script:
        def input_files = bam_files.collect { "--bam $it" }.join(' ')
        def exome_flag = (params.seq_type == 'WES') ? "--exome" : ""
    """
    /manta/bin/configManta.py \
        ${input_files} \
        --referenceFasta ${fasta} \
        ${exome_flag} \
        --runDir manta

    cd manta
    ./runWorkflow.py -j ${task.cpus}
    """
}