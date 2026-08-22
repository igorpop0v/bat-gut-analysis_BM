#!/usr/bin/env bash

# Export the ASV count table and GTDB taxonomy for analysis in R.

set -e

WORK_DIR="$1"
QIIME_IMAGE="quay.io/qiime2/qiime2:2026.7"
EXPORT_DIR="$WORK_DIR/exported_tables"

mkdir -p "$EXPORT_DIR"

echo "Exporting ASV count table..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:$WORK_DIR" \
  --workdir "$WORK_DIR" \
  "$QIIME_IMAGE" \
  qiime tools export \
    --input-path "$WORK_DIR/table_deblur_190.qza" \
    --output-path "$EXPORT_DIR/asv"

echo "Converting BIOM table to TSV..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:$WORK_DIR" \
  --workdir "$WORK_DIR" \
  "$QIIME_IMAGE" \
  biom convert \
    --input-fp "$EXPORT_DIR/asv/feature-table.biom" \
    --output-fp "$EXPORT_DIR/feature_table_asv.tsv" \
    --to-tsv

echo "Exporting GTDB taxonomy..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:$WORK_DIR" \
  --workdir "$WORK_DIR" \
  "$QIIME_IMAGE" \
  qiime tools export \
    --input-path "$WORK_DIR/taxonomy_gtdb_r232.qza" \
    --output-path "$EXPORT_DIR/taxonomy"

echo "Finished."
echo "ASV counts: $EXPORT_DIR/feature_table_asv.tsv"
echo "Taxonomy:   $EXPORT_DIR/taxonomy/taxonomy.tsv"
