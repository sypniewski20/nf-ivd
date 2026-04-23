# nf-ivd — Clinical Germline IVD Pipeline

A Nextflow DSL2 pipeline for **IVD-compliant germline variant calling** using GATK 4.6+ Best Practices in DRAGEN-mode. Achieves hardware-equivalent accuracy on standard CPU infrastructure via DRAGEN HMM and DRAGstr noise modelling.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Requirements](#requirements)
- [Repository Structure](#repository-structure)
- [Setup with Make](#setup-with-make)
- [Samplesheet Format](#samplesheet-format)
- [Quick Start](#quick-start)
- [Run Modes](#run-modes)
- [Configuration](#configuration)
- [Output Structure](#output-structure)
- [Benchmarking & Validation](#benchmarking--validation)
- [License & Compliance](#license--compliance)

---

## Overview

`nf-ivd` supports the following analysis types:

- **WGS & WES** sequencing inputs
- **Single-sample** germline variant calling
- **Joint/cohort calling** with pedigree-aware family priors (Trio support)
- **Structural variant (SV)** detection via Manta and gCNV
- **GIAB calibration mode** for benchmarking against NIST truth sets (HG002, HG003, HG004)

### TO DO

- **TUMOR and TUMOR/GERMLINE with Mutect2 + benchmark with SEQC2**
- **VEP annotation**
- **PacBio support? or seperate repo?**

---

## Pipeline Architecture

```
Input Layer
  ├── FASTQ mode  →  QC (FastQC + fastp)  →  Alignment (DragMap)
  │     ├── Standard local FASTQ
  │     └── NIST streaming mode (URL-based, MD5-tagged chunks)
  └── BAM mode    →  Skip to variant calling

Post-mapping Layer
  └── DRAGstr calibration  →  Coverage QC (Mosdepth)

Routing Layer (run_mode)
  ├── HC           →  HaplotypeCaller (DRAGEN-mode)  →  VQSR / filtering
  ├── SV           →  Manta  +  gCNV
  └── calibration  →  HC  →  hap.py benchmarking vs NIST truth

Reporting
  └── MultiQC  (FastQC + flagstat + Mosdepth + variant stats)
```

---

## Requirements

| Dependency | Version |
|---|---|
| Nextflow | ≥ 22.10 |
| Singularity / Apptainer | any recent |
| Container images | `core.sif`, `qc.sif`, `happi.sif`, `manta.sif` |
| Reference genome | GRCh38 FASTA + index |
| DRAGstr STR table | bundled with GATK 4.6+ |
| Pedigree file | optional (`.ped`, for joint-calling) |

> Container images must be version-locked in your local Singularity registry to satisfy IVD reproducibility requirements.

---

## Repository Structure

```
nf-ivd/
├── main.nf                         # Workflow entry point
├── nextflow.config                 # Parameters & profiles
├── fetch_ajtrio.R                  # Ashkenazim Trio GIAB data helper
├── modules/
│   └── functions.nf                # Samplesheet & BAM readers
├── subworkflows/
│   ├── mapping.nf                  # DragMap alignment
│   ├── qc.nf                       # FastQC, fastp, Mosdepth, NIST streaming QC
│   ├── multiqc.nf                  # Aggregated QC report
│   ├── HC.nf                       # HaplotypeCaller (DRAGEN-mode)
│   ├── manta.nf                    # Structural variant calling
│   ├── gCNV.nf                     # Copy number variant calling
│   └── calibration.nf              # GIAB benchmarking
└── deployment/
    ├── manifests/                  # Example samplesheets & PED files
    └── reference/truth/            # NIST v4.2.1 truth VCFs & BEDs
```

---

## Setup with Make

The `Makefile` automates the full environment setup — building containers, downloading references, and preparing interval lists. All artefacts land under `deployment/`.

### Targets

| Target | Description |
|---|---|
| `make all` | Runs `setup` + `data` (full environment bootstrap) |
| `make setup` | Runs `containers`, `strat`, `truth`, `fasta`, `add_resources` |
| `make containers` | Builds all four Singularity images from `.def` files (or pulls from Docker Hub for `happi.sif`) |
| `make strat` | Downloads & verifies GIAB stratification BEDs via `download_and_verify.R` |
| `make truth` | Downloads & verifies NIST v4.2.1 truth VCFs and BEDs |
| `make fasta` | Downloads the GRCh38 FASTA and builds the DragMap hash table |
| `make add_resources` | Fetches ploidy priors (GCS), 1000G PoN (GCS), Broad WGS intervals, ENCODE blacklist, and UCSC segmental duplications; then runs `refine_intervals.sh` |
| `make data` | Downloads Ashkenazim Trio FASTQ reads via `ashkenazim_trio_download.sh` (large — run manually) |
| `make clean` | Removes outputs, logs, and built `.sif` images |

### Container images built

| Image | Source |
|---|---|
| `core.sif` | `deployment/singularity/core.def` (GATK + DragMap) |
| `qc.sif` | `deployment/singularity/qc.def` (FastQC, fastp, Mosdepth, samtools) |
| `happi.sif` | `docker://mgibio/hap.py:v0.3.12` |
| `manta.sif` | `deployment/singularity/manta.def` |

### Typical first-time setup

```bash
# Build containers and download all references (requires gsutil + singularity)
make setup

# Download Ashkenazim Trio reads separately (large files — run when ready)
make data
```

> `make fasta` also runs `dragen-os --build-hash-table` via `core.sif`, so the FASTA directory must be writable and have sufficient disk space for the hash table.

---

## Samplesheet Format

The samplesheet is a comma-separated CSV passed to `--samplesheet`. The required columns differ by `--input_type`.

### FASTQ samplesheet (`--input_type fastq`)

Used by `Read_samplesheet()`. All columns are required unless marked optional.

| Column | Description | Required |
|---|---|---|
| `SM` | Sample name — used as the grouping key and BAM `SM` read group tag | ✅ |
| `ID` | Read group ID (`RGID`). Falls back to `{SM}_{LB}` if omitted | optional |
| `LB` | Library name (`RGLB`). Defaults to `unknown_lib` if omitted | optional |
| `PL` | Sequencing platform (`RGPL`), e.g. `ILLUMINA`. Defaults to `unknown_platform` | optional |
| `PU` | Platform unit (`RGPU`), e.g. flowcell barcode. Defaults to `unknown_unit` | optional |
| `R1` | Path or URL to Read 1 FASTQ (gzipped) | ✅ |
| `R2` | Path or URL to Read 2 FASTQ (gzipped) | ✅ |

Example:

```csv
SM,ID,LB,PL,PU,R1,R2
HG002,HG002_L001,lib1,ILLUMINA,HXXXXXX.1,/data/HG002_R1.fastq.gz,/data/HG002_R2.fastq.gz
HG003,HG003_L001,lib1,ILLUMINA,HXXXXXX.2,/data/HG003_R1.fastq.gz,/data/HG003_R2.fastq.gz
```

For **calibration mode**, `R1`/`R2` can be remote URLs (GIAB FTP/S3). The `ID` field is used as the MD5 chunk tag for streaming QC grouping.

### BAM samplesheet (`--input_type bam`)

Used by `Read_bam()`. All three columns are required.

| Column | Description | Required |
|---|---|---|
| `sampleID` | Sample identifier — used as the grouping key downstream | ✅ |
| `bam` | Path to sorted BAM file | ✅ |
| `bai` | Path to the corresponding BAM index (`.bai`) | ✅ |

Example:

```csv
sampleID,bam,bai
HG002,/data/HG002.bam,/data/HG002.bam.bai
HG003,/data/HG003.bam,/data/HG003.bam.bai
```

> Both samplesheet readers call `checkIfExists: true`, so the pipeline will fail fast at startup if any file path is invalid.

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/sypniewski20/nf-ivd.git
cd nf-ivd
```

### 2. Prepare `nextflow.config`

Set your local paths for the required parameters:

```groovy
params.fasta           = "/path/to/GRCh38.fa"
params.singularity_path = "/path/to/sif_images"
```

### 3. Run

**Single sample — WGS, FASTQ input:**

```bash
nextflow run main.nf \
    --input_type fastq \
    --run_mode HC \
    --samplesheet deployment/manifests/HG002_samplesheet.csv \
    -profile singularity
```

**Trio — joint calling with pedigree:**

```bash
nextflow run main.nf \
    --input_type fastq \
    --run_mode HC \
    --calling_mode cohort \
    --pedigree deployment/manifests/ashkenazim_trio.ped \
    --samplesheet deployment/manifests/trio_samplesheet.csv \
    -profile singularity
```

**Structural variants:**

```bash
nextflow run main.nf \
    --input_type bam \
    --run_mode SV \
    --samplesheet deployment/manifests/sample_bams.csv \
    -profile singularity
```

**GIAB calibration / benchmarking:**

```bash
nextflow run main.nf \
    --input_type fastq \
    --run_mode calibration \
    --samplesheet deployment/manifests/HG002_nist_samplesheet.csv \
    -profile singularity
```

---

## Run Modes

| `--run_mode` | Description |
|---|---|
| `HC` | HaplotypeCaller germline SNV/indel calling in DRAGEN-mode |
| `SV` | Structural variant calling (Manta + gCNV) |
| `calibration` | NIST streaming QC → HC → hap.py benchmarking vs GIAB truth |

Multiple modes can be combined as a comma-separated string, e.g. `--run_mode HC,SV`.

---

## Configuration

All parameters are set in `nextflow.config`. Key options:

| Parameter | Description | Default |
|---|---|---|
| `input_type` | Input data format: `fastq` or `bam` | `fastq` |
| `run_mode` | Analysis mode (see above) | `germline` |
| `seq_type` | `WGS` or `WES` — affects intervals & DRAGstr calibration | `WGS` |
| `calling_mode` | `single` or `cohort` (GenomicsDB / GenotypeGVCFs) | `cohort` |
| `dragen_mode` | Enable DRAGEN-equivalent HMM and parameters | `true` |
| `joint_calling` | Alias for `calling_mode cohort` | — |
| `pedigree` | Path to `.ped` file for Bayesian pedigree priors | `null` |
| `intervals_list` | Genomic regions for parallel scatter | `chr1..M` |
| `interval_padding` | Padding in bp around WES target intervals | `150` |
| `fasta` | Path to GRCh38 reference FASTA | `null` |
| `bed` | Target BED file (WES only) | `null` |
| `strat_dir` | GIAB stratification BEDs directory | `null` |
| `outfolder` | Root output directory | `results` |
| `singularity_path` | Path to directory containing `.sif` container images | `null` |

### Resource labels (Singularity profile)

| Label | CPUs | Memory | Wall time |
|---|---|---|---|
| `tiny` | 1 | 2 GB | 1 h |
| `small` | 2 | 6 GB | 4 h |
| `medium` | 4 | 16 GB | 12 h |
| `large` | 8 | 32 GB | 24 h |
| `xlarge` | 16 | 64 GB | 36 h |
| `xxlarge` | 32 | 128 GB | 48 h |

---

## Output Structure

Each run generates a timestamped directory under `--outfolder` (format: `YYYYMMDD_HHMMSS`):

```
results/
└── 20250423_120000/
    ├── bam/            # Sorted, indexed BAMs with full Read Group headers
    ├── vcf/            # Filtered final VCFs (SNV/indel and/or SV)
    ├── qc/             # Per-sample FastQC, fastp, flagstat, Mosdepth reports
    ├── multiqc/        # Aggregated MultiQC HTML report
    └── logs/
        ├── execution_timeline.html
        ├── execution_report.html
        ├── execution_trace.txt
        └── pipeline_dag.html
```

### Clinical integrity checks

- **BAM validation:** `samtools quickcheck` before any variant calling step.
- **MD5 checksums:** Generated for all final alignment files.
- **Variant normalisation:** All variants are decomposed, left-aligned, and normalised.
- **Full audit trail:** Nextflow timeline, trace, and DAG are always written to `logs/`.

---

## Benchmarking & Validation

The pipeline is pre-configured for benchmarking against the **Ashkenazim Trio (HG002 / HG003 / HG004)** using NIST v4.2.1 truth sets.

1. Download truth VCFs and BEDs (or use `fetch_ajtrio.R`) and place them under `deployment/reference/truth/`.
2. Point `strat_dir` at your local copy of the GIAB stratification BEDs.
3. Run with `--run_mode calibration`.
4. The pipeline produces `filtered_final.vcf.gz` ready for comparison with `hap.py` or `rtg-tools`.

Truth registry paths (configured in `nextflow.config`):

```groovy
truth {
    HG002 { vcf = "..."; bed = "..." }
    HG003 { vcf = "..."; bed = "..." }
    HG004 { vcf = "..."; bed = "..." }
}
```

---

## License & Compliance

Designed for research and clinical validation use. To satisfy IVD reproducibility requirements:

- Version-lock all container images (`core.sif`, `qc.sif`, `happi.sif`, `manta.sif`) in your local registry.
- Preserve the `logs/` directory for each run as your regulatory audit trail.
- Do not modify reference files between validation runs.
