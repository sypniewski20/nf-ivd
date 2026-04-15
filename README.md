# nf-ivd

Clinical-grade environment bootstrap for GIAB benchmarking and variant calling validation.

This repository **does not run the pipeline yet**.
It builds a **fully reproducible, auditable environment** in which the pipeline will run.

The goal is to guarantee:

* Reproducible containers
* Verified reference genome
* Audited data acquisition (reads and BEDs)
* Snapshot of environment state for traceability

---

## Repository Structure

```
.
├── Makefile
├── core.def
├── qc.def
├── happi.def
├── fasta_manifest.tsv
├── manifest.csv
├── giab_bed_manifest.csv
├── scripts/
│   ├── download_reads.R
│   ├── download_beds.R
│   ├── utils.R
│   └── init.R
```

---

## Entry Points (Important)

These are the only commands you need.

### 1. Full environment bootstrap

```
make all
```

This performs:

1. Build containers (core, qc, hap.py)
2. Download and verify reference genome
3. Initialize audit layer
4. Download reads and BED files
5. Create environment snapshot

---

### 2. Build only environment (no data)

```
make setup
```

Builds:

* `core.sif`
* `qc.sif`
* `happi.sif` (hap.py pinned by digest)
* Downloads and verifies reference from `fasta_manifest.tsv`
* Initializes audit ledger

---

### 3. Download all data (reads + beds)

```
make download
```

Uses:

* `manifest.csv` for reads
* `giab_bed_manifest.csv` for BED files

Each file is:

* Downloaded
* SHA256 verified
* Written to ledger

Safe to re-run. Already completed datasets are skipped.

---

### 4. Snapshot current environment state

```
make snapshot
```

Creates:

```
data/_run_snapshot.json
```

Contains:

* Timestamp
* Reference manifest used
* Containers used

This is the clinical traceability anchor.

---

## Reference Genome

Reference files are defined in:

```
fasta_manifest.tsv
```

Each entry contains:

```
file<TAB>sha256
```

Files are pulled from Broad GCS and verified before use.

---

## Containers

| Container | Source                      | Purpose                      |
| --------- | --------------------------- | ---------------------------- |
| core.sif  | core.def                    | aligners, samtools, bcftools |
| qc.sif    | qc.def                      | QC and metrics tools         |
| happi.sif | mgibio/hap.py pinned digest | GIAB benchmarking            |

`hap.py` is built from a **pinned Docker registry digest**, not a tag.

---

## Data Manifests

### Reads

```
manifest.csv
```

Columns:

* dataset
* accession (SRR/ERR) or
* source_url (http/ftp)

### BEDs

```
giab_bed_manifest.csv
```

Columns:

* dataset
* bed_type
* source_url

---

## Audit Layer

Every downloaded file gets:

* SHA256 computed
* Entry written into:

```
data/_ledger.json
```

This file is append-only and never overwritten.

---

## Cleaning

```
make clean
```

Removes:

* `fasta/`
* `data/`
* all `.sif` files
* inspect metadata
