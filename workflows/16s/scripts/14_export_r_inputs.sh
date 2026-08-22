#!/usr/bin/env bash

# Prepare clearly named 16S data files for downstream analysis in R.

set -e

WORK_DIR="$1"
QIIME_IMAGE="quay.io/qiime2/qiime2:2026.7"

EXPORT_DIR="$WORK_DIR/diversity/exports_for_r"
R_INPUT_DIR="$WORK_DIR/r_input"

mkdir -p "$EXPORT_DIR"
mkdir -p "$R_INPUT_DIR"

echo "Merging primary alpha-diversity metrics..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime metadata merge \
    --m-metadata1-file /data/diversity/core_metrics_18000/shannon_vector.qza \
    --m-metadata2-file /data/diversity/core_metrics_18000/simpson_vector.qza \
    --o-merged-metadata /data/diversity/alpha_shannon_simpson.qza

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime metadata merge \
    --m-metadata1-file /data/diversity/alpha_shannon_simpson.qza \
    --m-metadata2-file /data/diversity/core_metrics_18000/evenness_vector.qza \
    --o-merged-metadata /data/diversity/alpha_primary_metrics.qza

echo "Exporting alpha-diversity table..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime tools export \
    --input-path /data/diversity/alpha_primary_metrics.qza \
    --output-path /data/diversity/exports_for_r/alpha

cp \
  "$EXPORT_DIR/alpha/metadata.tsv" \
  "$R_INPUT_DIR/alpha_diversity_metrics_rarefied_18000.tsv"

# Replace the inconvenient QIIME 2 column names.
sed -i $'1c sample_id\tshannon\tsimpson\tpielou_evenness' \
  "$R_INPUT_DIR/alpha_diversity_metrics_rarefied_18000.tsv"

# Remove the optional QIIME 2 metadata type row.
sed -i '/^#q2:types/d' \
  "$R_INPUT_DIR/alpha_diversity_metrics_rarefied_18000.tsv"

echo "Exporting the unrarefied biological-sample ASV table..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime tools export \
    --input-path /data/diversity/table_biological.qza \
    --output-path /data/diversity/exports_for_r/asv_counts

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  biom convert \
    --input-fp /data/diversity/exports_for_r/asv_counts/feature-table.biom \
    --output-fp /data/r_input/asv_counts_unrarefied.tsv \
    --to-tsv

echo "Copying taxonomy artifact and metadata..."

cp \
  "$WORK_DIR/taxonomy_gtdb_r232.qza" \
  "$R_INPUT_DIR/asv_taxonomy_gtdb_r232.qza"

cp \
  "$WORK_DIR/diversity/metadata_biological.tsv" \
  "$R_INPUT_DIR/sample_metadata.tsv"

echo "Exporting the Bray-Curtis distance matrix..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime tools export \
    --input-path /data/diversity/core_metrics_18000/bray_curtis_distance_matrix.qza \
    --output-path /data/diversity/exports_for_r/bray_curtis_distance

cp \
  "$EXPORT_DIR/bray_curtis_distance/distance-matrix.tsv" \
  "$R_INPUT_DIR/bray_curtis_distance_matrix_rarefied_18000.tsv"

echo "Exporting Bray-Curtis PCoA results..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime tools export \
    --input-path /data/diversity/core_metrics_18000/bray_curtis_pcoa_results.qza \
    --output-path /data/diversity/exports_for_r/bray_curtis_pcoa

cp \
  "$EXPORT_DIR/bray_curtis_pcoa/ordination.txt" \
  "$R_INPUT_DIR/bray_curtis_pcoa_ordination_rarefied_18000.txt"

echo "Exporting PERMANOVA results..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime tools export \
    --input-path /data/diversity/bray_curtis_significance/bray_curtis_permanova.qzv \
    --output-path /data/diversity/exports_for_r/permanova

cp \
  "$EXPORT_DIR/permanova/permanova-pairwise.csv" \
  "$R_INPUT_DIR/bray_curtis_permanova_pairwise.csv"

echo "Exporting PERMDISP results..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime tools export \
    --input-path /data/diversity/bray_curtis_significance/bray_curtis_permdisp.qzv \
    --output-path /data/diversity/exports_for_r/permdisp

cp \
  "$EXPORT_DIR/permdisp/permdisp-pairwise.csv" \
  "$R_INPUT_DIR/bray_curtis_permdisp_pairwise.csv"

echo "Files prepared for R:"
ls -1 "$R_INPUT_DIR"
