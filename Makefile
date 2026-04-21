SINGULARITY ?= singularity
GSUTIL      ?= gsutil
WGET		?= wget
R           ?= Rscript

# ── Scripts ──────────────────────────────────────────────────────────────────────
RSCRIPT = Rscript

# ── Dirs ──────────────────────────────────────────────────────────────────────
REF_DIR    := deployment/reference
STRAT_DIR  := $(REF_DIR)/giab_stratifications
TRUTH_DIR  := $(REF_DIR)/truth
READS_DIR  := $(REF_DIR)/reads
FASTA_DIR  := $(REF_DIR)/fasta
ADD_RESOURCES    := $(REF_DIR)/additional_resources

# ── Manifests ─────────────────────────────────────────────────────────────────
STRAT_MANIFEST := deployment/manifests/giab_strat_manifest.csv
TRUTH_MANIFEST := deployment/manifests/truth_manifest.csv
READS_MANIFEST := deployment/manifests/reads_manifest.csv
FASTA_MANIFEST := deployment/manifests/fasta_manifest.csv

# ── Images ────────────────────────────────────────────────────────────────────

CORE_SIF  := deployment/singularity/core.sif
QC_SIF    := deployment/singularity/qc.sif
HAPPY_SIF := deployment/singularity/happi.sif
MANTA_SIF := deployment/singularity/manta.sif

# ── Additional Resources ───────────────────────────────────────────────────────
PLOIDY_PRIORS := gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.contig_ploidy_priors_homo_sapiens.tsv
PON_1K_GENOMES := gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
BROAD_INTERVALS := https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list
ENCODE_BLACKLIST := https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
UCSC_SEGDUPS := https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

########################################################

.PHONY: all setup containers \
        strat truth fasta \
        validate clean data

all: setup data

setup: containers strat truth fasta add_resources

# ── Containers ────────────────────────────────────────────────────────────────
containers: $(CORE_SIF) $(QC_SIF) $(HAPPY_SIF) $(MANTA_SIF)

$(CORE_SIF):
	$(SINGULARITY) build --fakeroot $@ deployment/singularity/core.def

$(QC_SIF):
	$(SINGULARITY) build --fakeroot $@ deployment/singularity/qc.def

$(HAPPY_SIF):
	$(SINGULARITY) build --disable-cache $@ docker://mgibio/hap.py:v0.3.12
 
$(MANTA_SIF):
	$(SINGULARITY) build --fakeroot $@ deployment/singularity/manta.def

# ── References ────────────────────────────────────────────────────────────────
strat:
	$(RSCRIPT) scripts/download_and_verify.R --manifest $(STRAT_MANIFEST) --dir $(STRAT_DIR) --snapshot_dir $(STRAT_DIR)

truth:
	$(RSCRIPT) scripts/download_and_verify.R --manifest $(TRUTH_MANIFEST) --dir $(TRUTH_DIR) --snapshot_dir $(TRUTH_DIR)

fasta:
	$(RSCRIPT) scripts/download_and_verify.R --manifest $(FASTA_MANIFEST) --dir $(FASTA_DIR) --snapshot_dir $(FASTA_DIR)

	# Build the hash table for the reference FASTA file
	$(SINGULARITY) run $(CORE_SIF) dragen-os --build-hash-table true \
											 --ht-reference ${FASTA_DIR}/*.fasta  \
											 --output-directory ${FASTA_DIR}

add_resources:
	$(GSUTIL) cp ${PLOIDY_PRIORS} ${ADD_RESOURCES}/
	$(GSUTIL) cp ${PON_1K_GENOMES} ${ADD_RESOURCES}/
	$(WGET) -P ${ADD_RESOURCES} ${BROAD_INTERVALS}
	$(WGET) -P ${ADD_RESOURCES} ${ENCODE_BLACKLIST}
	$(WGET) -P ${ADD_RESOURCES} ${UCSC_SEGDUPS}

	# Refine intervals

	$(SINGULARITY) run $(CORE_SIF) deployment/scripts/./refine_intervals.sh \
															$(ADD_RESOURCES)/wgs_calling_regions.hg38.interval_list \
															$(ADD_RESOURCES)/ENCFF356LFX.bed.gz \
															$(ADD_RESOURCES)/genomicSuperDups.txt.gz \
															$(FASTA_DIR)/*.dict \
															$(ADD_RESOURCES)


# ── Reads (manual step — too large for default pipeline) ─────────────────────
data:
	scripts/ashkenazim_trio_download.sh ${READS_DIR}

# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	rm -rf $(OUTPUT_DIR) reference logs singularity/*.sif