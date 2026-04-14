# nf-ivd

Reproducible genomics environment for clinical-style benchmarking workflows using:

- DRAGMAP-compatible GRCh38 reference
- GIAB / SEQC2 / COLO829 benchmark datasets
- SHA256-verified ingestion layer
- Singularity containerized tools
- Run snapshot audit layer

---

# 🧬 Execution Model

The system is structured into three deterministic layers:

```

ENVIRONMENT → DOWNLOAD → SNAPSHOT

````

- **ENVIRONMENT**: builds tools + reference
- **DOWNLOAD**: acquires datasets with integrity checks
- **SNAPSHOT**: captures immutable run state for auditability

---

# 🚀 Full Entry Points (IMPORTANT)

## 1. Full system bootstrap

```bash
make all
````

Runs full reproducible setup:

* container build (core, qc, hap.py)
* reference download + SHA256 verification
* dataset download (reads + BEDs)
* run snapshot generation

---

## 2. Environment setup only

```bash
make setup
```

Includes:

* Singularity container builds
* DRAGMAP-compatible reference download
* reference SHA256 verification
* audit directory initialization

---

## 3. Data download only

```bash
make download
```

Includes:

* SRA / FASTQ acquisition
* direct HTTP/FTP dataset download
* GIAB BED truth set download
* SHA256 validation per file
* resume-safe execution (ledger-based)

---

## 4. Snapshot only

```bash
make snapshot
```

Generates:

```
data/_run_snapshot.json
```

Contains:

* timestamp
* container identifiers
* reference linkage
* execution state snapshot

---

# 📦 Data Model

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

# 🧬 Reference Genome

Based on:

* GRCh38 (Broad / DRAGMAP compatible)
* pinned via `reference_manifest.tsv`
* SHA256-verified downloads from GCS

Core files:

* FASTA
* FAI index
* dictionary (.dict)
* STR metadata (.str)

---

# 📁 Output Structure

```
reference/         # immutable reference assets
data/              # downloaded datasets
data/_ledger.json  # download tracking
data/_run_snapshot.json  # run audit snapshot
```

---

# 📟 Command Reference

| Command         | Purpose                 |
| --------------- | ----------------------- |
| `make all`      | Full system bootstrap   |
| `make setup`    | Build environment only  |
| `make download` | Acquire datasets        |
| `make snapshot` | Generate audit snapshot |
| `make clean`    | Reset workspace         |

```
