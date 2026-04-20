## ⚙️ Infrastructure & Setup (Makefile)

The included `Makefile` manages the clinical environment's lifecycle, ensuring that all reference data and containers are cryptographically verified before use.

### 1. Environment Initialization
To build the required clinical containers and fetch reference data:

```bash
make setup
```

### 2. Available Commands

| Command | Description |
| :--- | :--- |
| `make containers` | Builds Singularity images (`core.sif`, `qc.sif`) from local `.def` files. |
| `make fasta` | Downloads and verifies the GRCh38 reference genome via manifest. |
| `make truth` | Fetches GIAB/NIST truth sets for HG002-004 and SEQC2 benchmarking. |
| `make strat` | Downloads GIAB genomic stratifications for specific region analysis. |
| `make pon` | Synchronizes the GATK Best Practices Panel of Normals (PoN). |
| `make data` | Triggers the high-volume benchmark FASTQ download script. |
| `make clean` | Wipes the `output/`, `reference/`, and `logs/` to reset the environment. |

### 3. Data Integrity & Verification
To maintain **IVD (In Vitro Diagnostics)** compliance, all data downloads are governed by manifests located in the `manifests/` directory. 

The `setup` process calls `scripts/download_and_verify.R`, which:
* Compares downloaded file hashes against the `manifests/*.csv`.
* Creates a `snapshot_dir` to ensure that reference data cannot be modified post-download without triggering a validation error.

### 4. Container Security
The `core.sif` and `qc.sif` images are built using the `--fakeroot` flag from local definitions to ensure no unverified layers are introduced from external registries.