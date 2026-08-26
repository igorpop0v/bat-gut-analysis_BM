#!/usr/bin/env bash

# Predict gene families and metabolic pathways from 16S ASV data.

set -e

PICRUST_DIR="$1"
THREADS="${2:-12}"

PICRUST_IMAGE="quay.io/biocontainers/picrust2:2.6.3--pyhdfd78af_2"
OUTPUT_DIR="$PICRUST_DIR/output/picrust2_v2.6.3"

mkdir -p "$PICRUST_DIR/output"

echo "Running PICRUSt2 with $THREADS threads..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$PICRUST_DIR:/picrust2" \
  --workdir /picrust2 \
  "$PICRUST_IMAGE" \
  picrust2_pipeline.py \
    --study_fasta /picrust2/input/rep_seqs/dna-sequences.fasta \
    --input /picrust2/input/asv_table/feature-table.biom \
    --output /picrust2/output/picrust2_v2.6.3 \
    --processes "$THREADS" \
    --verbose

echo "PICRUSt2 analysis completed."
echo "Results: $OUTPUT_DIR"
