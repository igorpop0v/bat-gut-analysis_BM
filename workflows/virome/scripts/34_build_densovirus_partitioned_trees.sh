#!/usr/bin/env bash

set -euo pipefail

INPUT_DIR="$1"
OUTPUT_DIR="$2"
THREADS="$3"

SELECTED_DATASET="${4:-all}"

IQTREE_IMAGE="evolbioinfo/iqtree:v3.1.1"
RANDOM_SEED="20260831"

if [[ "${SELECTED_DATASET}" == "all" ]]
then
    DATASETS=(
        "primary"
        "study_inclusive"
        "sensitivity"
    )
else
    DATASETS=(
        "${SELECTED_DATASET}"
    )
fi

mkdir -p "${OUTPUT_DIR}"

for DATASET in "${DATASETS[@]}"
do
    echo
    echo "Building partitioned NSP1 + SP1 tree: ${DATASET}..."

    mkdir -p "${OUTPUT_DIR}/${DATASET}"

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      --volume "${INPUT_DIR}:/input:ro" \
      --volume "${OUTPUT_DIR}:/output" \
      "${IQTREE_IMAGE}" \
      -s "/input/${DATASET}/densovirus_${DATASET}_nsp1_sp1.faa" \
      -p "/input/${DATASET}/densovirus_${DATASET}_partitions.nex" \
      -st AA \
      -m MFP \
      -B 1000 \
      -alrt 1000 \
      -bnni \
      -T "${THREADS}" \
      -seed "${RANDOM_SEED}" \
      --prefix \
      "/output/${DATASET}/densovirus_${DATASET}_nsp1_sp1"
done

echo
echo "Partitioned IQ-TREE analyses completed."
echo "Results: ${OUTPUT_DIR}"
