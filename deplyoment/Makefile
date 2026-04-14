SINGULARITY ?= singularity
GSUTIL ?= gsutil
R ?= Rscript

########################################################
# IMAGES
########################################################

CORE_DEF = core.def
CORE_SIF = core.sif

QC_DEF = qc.def
QC_SIF = qc.sif

HAPPY_DEF = happi.def
HAPPY_SIF = happi.sif

########################################################
# REFERENCE
########################################################

REF_DIR = reference
REF_MANIFEST = reference_manifest.tsv
REF_LEDGER = $(REF_DIR)/_ledger.json

GCS_BASE = gs://gcp-public-data--broad-references/hg38/v0

########################################################
# DATA
########################################################

DATA_DIR = data

READS_SCRIPT = scripts/download_reads.R
BEDS_SCRIPT  = scripts/download_beds.R

READS_LEDGER = $(DATA_DIR)/_ledger.json
BEDS_LEDGER  = $(DATA_DIR)/_ledger.json
RUN_SNAPSHOT = $(DATA_DIR)/_run_snapshot.json

########################################################
# DEFAULT
########################################################

all: setup download snapshot

########################################################
# ======================================================
# FULL SYSTEM SETUP (NEW IMPORTANT ENTRYPOINT)
# ======================================================
########################################################

setup: env

env: containers reference audit

containers: core qc happi

core:
	@echo "[ENV] Building CORE container"
	$(SINGULARITY) build --fakeroot $(CORE_SIF) $(CORE_DEF)

qc:
	@echo "[ENV] Building QC container"
	$(SINGULARITY) build --fakeroot $(QC_SIF) $(QC_DEF)

happi:
	@echo "[ENV] Building hap.py container"
	$(SINGULARITY) build --fakeroot $(HAPPY_SIF) $(HAPPY_DEF)

reference: ref-download ref-verify ref-ledger

ref-download:
	@echo "[ENV] Downloading reference"
	mkdir -p $(REF_DIR)

	@tail -n +2 $(REF_MANIFEST) | while read f sha; do \
		echo "[REF] $$f"; \
		$(GSUTIL) cp $(GCS_BASE)/$$f $(REF_DIR)/$$f; \
	done

ref-verify:
	@echo "[ENV] Verifying reference integrity"
	cd $(REF_DIR) && \
	while read f sha; do \
		[ "$$f" = "file" ] && continue; \
		echo "$$sha  $$f" | sha256sum -c - || exit 1; \
	done < ../$(REF_MANIFEST)

ref-ledger:
	@test -f $(REF_LEDGER) || echo "[]" > $(REF_LEDGER)

audit:
	@echo "[ENV] Initializing audit layer"
	mkdir -p $(DATA_DIR)
	@test -f $(READS_LEDGER) || echo "[]" > $(READS_LEDGER)

########################################################
# ======================================================
# DOWNLOAD LAYER
# ======================================================
########################################################

download: download-all

download-all: reads beds

reads:
	@echo "[DOWNLOAD] Reads acquisition"
	mkdir -p $(DATA_DIR)
	$(R) $(READS_SCRIPT)

beds:
	@echo "[DOWNLOAD] BED acquisition"
	mkdir -p $(DATA_DIR)
	$(R) $(BEDS_SCRIPT)

########################################################
# ======================================================
# SNAPSHOT LAYER (CLINICAL AUDIT CORE)
# ======================================================
########################################################

snapshot:
	@echo "[SNAPSHOT] Generating run snapshot"
	mkdir -p $(DATA_DIR)

	@echo "{" > $(RUN_SNAPSHOT)
	@echo "  \"timestamp\": \"$(shell date -Iseconds)\"," >> $(RUN_SNAPSHOT)
	@echo "  \"reference_manifest\": \"$(REF_MANIFEST)\"," >> $(RUN_SNAPSHOT)
	@echo "  \"containers\": {" >> $(RUN_SNAPSHOT)
	@echo "    \"core\": \"$(CORE_SIF)\"," >> $(RUN_SNAPSHOT)
	@echo "    \"qc\": \"$(QC_SIF)\"," >> $(RUN_SNAPSHOT)
	@echo "    \"hap.py\": \"$(HAPPY_SIF)\"" >> $(RUN_SNAPSHOT)
	@echo "  }" >> $(RUN_SNAPSHOT)
	@echo "}" >> $(RUN_SNAPSHOT)

########################################################
# CLEAN
########################################################

clean:
	rm -rf $(REF_DIR) $(DATA_DIR) *.sif