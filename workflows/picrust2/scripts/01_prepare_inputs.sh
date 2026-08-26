#!/usr/bin/env bash

# Prepare representative sequences and the unrarefied ASV table
# from biological samples for PICRUSt2.

set -e

QIIME_DIR="$1"
PICRUST_DIR="$2"

QIIME_IMAGE="quay.io/qiime2/qiime2:2026.7"

mkdir -p "$PICRUST_DIR/input"

echo "Filtering representative sequences..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$QIIME_DIR:/qiime" \
  --volume "$PICRUST_DIR:/picrust2" \
  "$QIIME_IMAGE" \
  qiime feature-table filter-seqs \
    --i-data /qiime/rep_seqs_deblur_190.qza \
    --i-table /qiime/diversity/table_biological.qza \
    --o-filtered-data /picrust2/input/rep_seqs_biological.qza

echo "Exporting representative sequences..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$PICRUST_DIR:/picrust2" \
  "$QIIME_IMAGE" \
  qiime tools export \
    --input-path /picrust2/input/rep_seqs_biological.qza \
    --output-path /picrust2/input/rep_seqs

echo "Exporting the unrarefied biological-sample ASV table..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$QIIME_DIR:/qiime" \
  --volume "$PICRUST_DIR:/picrust2" \
  "$QIIME_IMAGE" \
  qiime tools export \
    --input-path /qiime/diversity/table_biological.qza \
    --output-path /picrust2/input/asv_table

echo "PICRUSt2 input files:"
ls -lh \
  "$PICRUST_DIR/input/rep_seqs/dna-sequences.fasta" \
  "$PICRUST_DIR/input/asv_table/feature-table.biom"
