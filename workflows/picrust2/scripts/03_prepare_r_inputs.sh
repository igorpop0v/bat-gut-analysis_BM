#!/usr/bin/env bash

# Add MetaCyc pathway descriptions and prepare PICRUSt2 files for R.

set -e

PICRUST_DIR="$1"
QIIME_DIR="$2"

PICRUST_IMAGE="quay.io/biocontainers/picrust2:2.6.3--pyhdfd78af_2"

RESULT_DIR="$PICRUST_DIR/output/picrust2_v2.6.3"
R_INPUT_DIR="$PICRUST_DIR/r_input"

mkdir -p "$R_INPUT_DIR"

echo "Adding MetaCyc pathway descriptions..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$PICRUST_DIR:/picrust2" \
  "$PICRUST_IMAGE" \
  add_descriptions.py \
    --input /picrust2/output/picrust2_v2.6.3/pathways_out/path_abun_unstrat.tsv.gz \
    --map_type METACYC \
    --output /picrust2/r_input/metacyc_pathway_abundance_with_descriptions.tsv.gz

echo "Copying PICRUSt2 QC files..."

cp \
  "$RESULT_DIR/pathways_out/path_abun_unstrat.tsv.gz" \
  "$R_INPUT_DIR/metacyc_pathway_abundance.tsv.gz"

cp \
  "$RESULT_DIR/combined_marker_predicted_and_nsti.tsv.gz" \
  "$R_INPUT_DIR/asv_nsti.tsv.gz"

cp \
  "$RESULT_DIR/EC_metagenome_out/weighted_nsti.tsv.gz" \
  "$R_INPUT_DIR/sample_weighted_nsti.tsv.gz"

echo "Copying metadata and original ASV counts..."

cp \
  "$QIIME_DIR/diversity/metadata_biological.tsv" \
  "$R_INPUT_DIR/sample_metadata.tsv"

cp \
  "$QIIME_DIR/r_input/asv_counts_unrarefied.tsv" \
  "$R_INPUT_DIR/asv_counts_unrarefied.tsv"

echo "Files prepared for R:"
ls -lh "$R_INPUT_DIR"
