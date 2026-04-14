#!/usr/bin/env bash
set -euo pipefail

MANIFEST="build_manifest.tsv"

# ---- optional: validate manifest integrity ----
# echo "EXPECTED_SHA  $MANIFEST" | sha256sum -c -

declare -A VERSIONS

while IFS=$'\t' read -r tool version url sha; do
    VERSIONS["$tool"]="$version"
done < "$MANIFEST"

require() {
    if [[ -z "${VERSIONS[$1]:-}" ]]; then
        echo "ERROR: missing version for $1" >&2
        exit 1
    fi
}

require GATK
require DRAGMAP
require SAMTOOLS
require HAPY

TAG="pon-benchmark-gatk-${VERSIONS[GATK]}"

docker build \
  --build-arg GATK_VERSION="${VERSIONS[GATK]}" \
  --build-arg DRAGMAP_VERSION="${VERSIONS[DRAGMAP]}" \
  --build-arg SAMTOOLS_VERSION="${VERSIONS[SAMTOOLS]}" \
  --build-arg HAPY_VERSION="${VERSIONS[HAPY]}" \
  -t "$TAG" .