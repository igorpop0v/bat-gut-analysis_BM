#!/usr/bin/env bash

set -euo pipefail

REFERENCE_DIR="$1"
STUDY_DIR="$2"
OUTPUT_DIR="$3"
THREADS="$4"

MAFFT_IMAGE="evolbioinfo/mafft:v7.526"

mkdir -p \
  "${OUTPUT_DIR}/combined" \
  "${OUTPUT_DIR}/mafft"

echo "Combining primary NSP1 sequences..."

awk '1' \
  "${REFERENCE_DIR}/external_nsp1_primary.faa" \
  "${STUDY_DIR}/study_nsp1_primary.faa" \
  > "${OUTPUT_DIR}/combined/nsp1_primary.faa"

echo "Combining sensitivity NSP1 sequences..."

awk '1' \
  "${REFERENCE_DIR}/external_nsp1_sensitivity.faa" \
  "${STUDY_DIR}/study_nsp1_sensitivity.faa" \
  > "${OUTPUT_DIR}/combined/nsp1_sensitivity.faa"

echo "Combining primary SP1 sequences..."

awk '1' \
  "${REFERENCE_DIR}/external_sp1_primary.faa" \
  "${STUDY_DIR}/study_sp1_primary.faa" \
  > "${OUTPUT_DIR}/combined/sp1_primary.faa"

echo "Combining sensitivity SP1 sequences..."

awk '1' \
  "${REFERENCE_DIR}/external_sp1_sensitivity.faa" \
  "${STUDY_DIR}/study_sp1_sensitivity.faa" \
  > "${OUTPUT_DIR}/combined/sp1_sensitivity.faa"

run_mafft() {

    DATASET="$1"

    echo "Aligning ${DATASET}..."

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      --volume "${OUTPUT_DIR}:/data" \
      "${MAFFT_IMAGE}" \
      --localpair \
      --maxiterate 1000 \
      --thread "${THREADS}" \
      "/data/combined/${DATASET}.faa" \
      > "${OUTPUT_DIR}/mafft/${DATASET}_mafft.faa"
}

run_mafft "nsp1_primary"
run_mafft "nsp1_sensitivity"
run_mafft "sp1_primary"
run_mafft "sp1_sensitivity"

echo
echo "Sequence counts in the resulting alignments:"

for ALIGNMENT in "${OUTPUT_DIR}"/mafft/*.faa
do
    printf "%s\t" "$(basename "${ALIGNMENT}")"
    grep -c '^>' "${ALIGNMENT}"
done

echo
echo "MAFFT alignments completed."
echo "Results: ${OUTPUT_DIR}/mafft"
