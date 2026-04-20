# Clinical Germline IVD Pipeline (GATK-DRAGEN Mode)

## Overview
This pipeline is a high-fidelity genomic workflow designed for **IVD-compliant** germline variant calling. It implements the **GATK 4.6+ Best Practices** using the **DRAGEN-mode** HMM and **DRAGstr** noise modeling to achieve hardware-equivalent accuracy on standard CPU infrastructure.

The workflow supports:
* **WGS & WES** sequencing types.
* **Single-Sample** analysis for routine diagnostics.
* **Joint-Calling (Trio/Cohort)** for pedigree-aware family analysis.
* **Clinical Benchmarking** via GIAB (Genome in a Bottle) truth sets.

---

## Architecture
The pipeline is built with modularity as a core principle:
* **Input Layer:** Handles FASTQ-to-BAM or BAM-to-VCF (checkpoint) entry points.
* **DRAGstr Calibration:** Mandatory per-sample noise model calibration.
* **Variant Calling:** Parallelized scattering across genomic intervals.
* **Refinement:** Bayesian pedigree refinement using `.ped` files.

---

## Quick Start

### 1. Requirements
* Nextflow (>= 22.10)
* Singularity / Apptainer
* References: GRCh38 Fasta, DRAGstr STR Table, and Pedigree file (optional).

### 2. Execution

**Run as a Single Sample (WGS):**
```bash
nextflow run main.nf \
    --input_type fastq \
    --run_mode HC \
    --joint_calling false \
    --samplesheet deployment/manifests/HG002_samplesheet.csv \
    -profile singularity
```

**Run as a Trio (Joint Calling):**
```bash
nextflow run main.nf \
    --input_type fastq \
    --run_mode HC \
    --joint_calling true \
    --pedigree deployment/manifests/ashkenazim_trio.ped \
    --samplesheet deployment/manifests/trio_samplesheet.csv \
    -profile singularity
```

## Configuration
Key parameters in `nextflow.config`:

| Parameter | Description | Default |
| :--- | :--- | :--- |
| `joint_calling` | Enable GenomicsDB & GenotypeGVCFs branch | `true` |
| `seq_type` | `WGS` or `WES` (affects intervals & calibration) | `WGS` |
| `dragen_mode` | Enforces DRAGEN-equivalent HMM and parameters | `true` |
| `intervals_list` | List of genomic regions for parallel scattering | `chr1..M` |
| `pedigree` | Path to the validated .ped file for family priors | `null` |

---

## Clinical Traceability & Data Integrity
Every run generates a timestamped `runID` folder (Format: `YYYYMMDD_HHMMSS`) containing:

1.  **VCFs:** Filtered and (optionally) Pedigree-refined variant calls.
2.  **BAMs:** Sorted, indexed, and MD5-verified alignments with full Read Group headers.
3.  **QC Reports:** Integrated MultiQC report (FastQC, Flagstat, Mosdepth).
4.  **Audit Logs:** Nextflow timeline, trace, and execution DAG for regulatory compliance.

### Integrity Checks
* **BAM Validation:** `samtools quickcheck` is performed on all alignments before variant calling.
* **File Traceability:** MD5 checksums are generated for all final alignment files.
* **Normalization:** All variants are decomposed and normalized (left-aligned) for clinical consistency.

---

## Validation & Benchmarking
This workflow is pre-configured to benchmark against the **Ashkenazim Trio (HG002, HG003, HG004)**. 

To run a validation study:
1. Ensure the `truth` registry in `nextflow.config` points to your local copies of the NIST v4.2.1 benchmarks.
2. The pipeline will produce a `filtered_final.vcf.gz` which can be compared using `hap.py` or `rtg-tools`.

---

## License & Compliance
Designed for research and clinical validation. Ensure all container images (`core.sif`, `qc.sif`) are version-locked in your local registry to maintain IVD reproducibility.