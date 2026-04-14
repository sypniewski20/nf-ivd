# nf-ivd

Reproducible genomics environment for benchmarking workflows using:

- DRAGMAP-compatible GRCh38 reference
- GIAB / SEQC2 / COLO829 benchmark datasets
- SHA256-verified ingestion layer
- Singularity containerized tools
- Run snapshot audit layer

---

# Execution Model

The system is structured into three deterministic layers:

```

ENVIRONMENT → DOWNLOAD → SNAPSHOT

````

- ENVIRONMENT: builds tools and reference assets
- DOWNLOAD: acquires datasets with integrity checks
- SNAPSHOT: captures immutable execution state for auditability

---

# Full Entry Points

## 1. Full system bootstrap

```bash
make all
````

Runs:

* container build (core, qc, hap.py)
* reference download + SHA256 verification
* dataset download (reads + BEDs)
* run snapshot generation

---

## 2. Environment setup only

```bash id="setup_cmd"
make setup
```

Includes:

* Singularity container builds
* DRAGMAP-compatible reference download
* reference verification using `reference_manifest.tsv`
* audit directory initialization

---

## 3. Data download only

```bash id="download_cmd"
make download
```

Includes:

* SRA / FASTQ acquisition
* direct HTTP/FTP dataset download
* GIAB BED truth set download
* SHA256 validation per file
* resume-safe execution using ledger files

---

## 4. Snapshot only

```bash id="snapshot_cmd"
make snapshot
```

Generates:

```
data/_run_snapshot.json
```

Contains:

* timestamp
* container identifiers
* reference manifest linkage
* execution state snapshot

---

# Data Model

## Reads

Supported:

* SRA (SRR / ERR)
* direct FASTQ download

Datasets:

* GIAB NA12878 / HG002
* SEQC2 tumor/normal
* COLO829 melanoma benchmark

## BEDs

Includes:

* callable regions
* excluded regions
* stratification BEDs (GIAB)

---

# Reference Genome

Based on:

* GRCh38 (Broad / DRAGMAP compatible)
* pinned via `reference_manifest.tsv`
* SHA256 verification against manifest

Core files:

* Homo_sapiens_assembly38.fasta
* Homo_sapiens_assembly38.fasta.fai
* Homo_sapiens_assembly38.dict
* Homo_sapiens_assembly38.str

---

# Output Structure

```
reference/                  # immutable reference assets
data/                       # downloaded datasets
reference_manifest.tsv     # reference file registry
giab_bed_manifest.csv      # BED registry
manifest.csv              # read dataset registry

data/_ledger.json         # download tracking
data/_run_snapshot.json   # run audit snapshot
```

---

# Command Reference

| Command       | Purpose                 |
| ------------- | ----------------------- |
| make all      | Full system bootstrap   |
| make setup    | Build environment only  |
| make download | Acquire datasets        |
| make snapshot | Generate audit snapshot |
| make clean    | Reset workspace         |

```

