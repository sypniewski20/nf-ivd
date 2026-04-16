READS_DIR=$1

# Ensure the directory exists first
mkdir -p ${READS_DIR}

java -jar /home/mateuszsypniewski/ENA_downloader/ena-file-downloader.jar \
--accessions=SRR14724532,SRR14724531,SRR14724530,SRR7890827,SRR7890824 \
--format=READS_FASTQ \
--protocol=ftp \
--destination=${READS_DIR}