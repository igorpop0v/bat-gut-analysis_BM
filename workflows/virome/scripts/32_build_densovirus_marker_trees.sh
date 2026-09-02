#!/usr/bin/env bash

set -euo pipefail

ALIGNMENT_DIR="$1"
OUTPUT_DIR="$2"
THREADS="$3"

IQTREE_IMAGE="evolbioinfo/iqtree:v3.1.1"
RANDOM_SEED="20260831"

DATASETS=(
    "nsp1_primary"
    "nsp1_sensitivity"
    "sp1_primary"
    "sp1_sensitivity"
)

mkdir -p "${OUTPUT_DIR}"

for DATASET in "${DATASETS[@]}"
do
    echo
    echo "Building Maximum Likelihood tree for ${DATASET}..."

    mkdir -p "${OUTPUT_DIR}/${DATASET}"

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      --volume "${ALIGNMENT_DIR}:/input:ro" \
      --volume "${OUTPUT_DIR}:/output" \
      "${IQTREE_IMAGE}" \
      -s "/input/${DATASET}_gt70.faa" \
      -st AA \
      -m MFP \
      -B 1000 \
      -alrt 1000 \
      -bnni \
      -T "${THREADS}" \
      -seed "${RANDOM_SEED}" \
      --prefix "/output/${DATASET}/${DATASET}"
done

echo
echo "IQ-TREE analyses completed."
echo "Results: ${OUTPUT_DIR}"
