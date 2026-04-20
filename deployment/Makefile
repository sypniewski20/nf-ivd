SINGULARITY ?= singularity
GSUTIL      ?= gsutil
R           ?= Rscript

# ── Scripts ──────────────────────────────────────────────────────────────────────
RSCRIPT = Rscript

# ── Dirs ──────────────────────────────────────────────────────────────────────
STRAT_DIR  := reference/giab_stratifications
TRUTH_DIR  := reference/truth
READS_DIR  := reference/reads
FASTA_DIR  := reference/fasta
REF_DIR    := reference
PON_DIR    := reference/pon
OUTPUT_DIR := output
LOG_DIR    := logs

# ── Manifests ─────────────────────────────────────────────────────────────────
STRAT_MANIFEST := manifests/giab_strat_manifest.csv
TRUTH_MANIFEST := manifests/truth_manifest.csv
READS_MANIFEST := manifests/reads_manifest.csv
FASTA_MANIFEST := manifests/fasta_manifest.csv

# ── Images ────────────────────────────────────────────────────────────────────

CORE_SIF  := singularity/core.sif
QC_SIF    := singularity/qc.sif
HAPPY_SIF := singularity/happi.sif
MANTA_SIF := singularity/manta.sif

########################################################

.PHONY: all setup containers \
        strat truth fasta \
        validate clean data

all: setup data

setup: containers strat truth fasta pon

# ── Containers ────────────────────────────────────────────────────────────────
containers: $(CORE_SIF) $(QC_SIF) $(HAPPY_SIF) $(MANTA_SIF)

$(CORE_SIF):
	$(SINGULARITY) build --fakeroot $@ singularity/core.def

$(QC_SIF):
	$(SINGULARITY) build --fakeroot $@ singularity/qc.def

$(HAPPY_SIF):
	$(SINGULARITY) build --disable-cache $@ docker://mgibio/hap.py:v0.3.12
 
$(MANTA_SIF):
	$(SINGULARITY) build --fakeroot $@ singularity/manta.def

# ── References ────────────────────────────────────────────────────────────────
strat:
	$(RSCRIPT) scripts/download_and_verify.R --manifest $(STRAT_MANIFEST) --dir $(STRAT_DIR) --snapshot_dir $(STRAT_DIR)

truth:
	$(RSCRIPT) scripts/download_and_verify.R --manifest $(TRUTH_MANIFEST) --dir $(TRUTH_DIR) --snapshot_dir $(TRUTH_DIR)

fasta:
	$(RSCRIPT) scripts/download_and_verify.R --manifest $(FASTA_MANIFEST) --dir $(FASTA_DIR) --snapshot_dir $(FASTA_DIR)

	# Build the hash table for the reference FASTA file
	${SINGULARITY} run $(CORE_SIF) dragen-os --build-hash-table true \
											 --ht-reference ${FASTA_DIR}/*.fasta  \
											 --output-directory ${FASTA_DIR}

pon:
#TODO: add PON download and verification
# 	gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz $(PON_DIR)/.

# ── Reads (manual step — too large for default pipeline) ─────────────────────
data:
	scripts/ashkenazim_trio_download.sh ${READS_DIR}

# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	rm -rf $(OUTPUT_DIR) reference logs singularity/*.sif