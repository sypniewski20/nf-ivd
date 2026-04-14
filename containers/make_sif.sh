#!/usr/bin/env bash
set -euo pipefail

GATK_VERSION=$(awk -F'\t' '$1=="GATK"{print $2}' build_manifest.tsv)
IMAGE="pon-benchmark-gatk${GATK_VERSION}"

singularity build ${IMAGE}.sif docker-daemon://${IMAGE}