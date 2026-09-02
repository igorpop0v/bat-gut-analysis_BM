#!/usr/bin/env bash

set -euo pipefail

INPUT_DIR="$1"
OUTPUT_DIR="$2"

TRIMAL_IMAGE="evolbioinfo/trimal:v1.5.1"
MINIMUM_OCCUPANCY="0.7"

mkdir -p "${OUTPUT_DIR}"

DATASETS=(
    "nsp1_primary"
    "nsp1_sensitivity"
    "sp1_primary"
    "sp1_sensitivity"
)

for DATASET in "${DATASETS[@]}"
do
    echo "Trimming ${DATASET} with minimum occupancy ${MINIMUM_OCCUPANCY}..."

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      --volume "${INPUT_DIR}:/input:ro" \
      --volume "${OUTPUT_DIR}:/output" \
      "${TRIMAL_IMAGE}" \
      -in "/input/${DATASET}_mafft.faa" \
      -out "/output/${DATASET}_gt70.faa" \
      -gt "${MINIMUM_OCCUPANCY}" \
      -fasta
done

echo
echo "Trimmed sequence counts:"

for ALIGNMENT in "${OUTPUT_DIR}"/*.faa
do
    printf "%s\t" "$(basename "${ALIGNMENT}")"
    grep -c '^>' "${ALIGNMENT}"
done

echo
echo "Fixed-threshold trimming completed."
echo "Results: ${OUTPUT_DIR}"
