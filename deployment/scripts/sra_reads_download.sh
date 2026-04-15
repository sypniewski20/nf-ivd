READS_MANIFEST=$1
SRA=$( tail -n +2 ${READS_MANIFEST} | cut -d, -f3 )

mkdir -p reference/reads

for sra in $SRA; do
    
    mkdir -p reference/reads/${sra}
    prefetch --max-size 100G -p ${sra} --output-directory reference/reads/${sra}
    #fastq-dump --gzip --split-files --outdir reference/reads/${sra} reference/reads/${sra}/${sra}.sra

done