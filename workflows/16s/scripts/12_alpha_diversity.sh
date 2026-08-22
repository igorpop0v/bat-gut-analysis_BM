#!/usr/bin/env bash

# Calculate Simpson diversity and create alpha-diversity
# group-significance visualizations.

set -e

WORK_DIR="$1"
QIIME_IMAGE="quay.io/qiime2/qiime2:2026.7"
CORE_DIR="$WORK_DIR/diversity/core_metrics_18000"
OUTPUT_DIR="$WORK_DIR/diversity/alpha_significance"

mkdir -p "$OUTPUT_DIR"

echo "Calculating Simpson diversity..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime diversity alpha \
    --i-table /data/diversity/core_metrics_18000/rarefied_table.qza \
    --p-metric simpson \
    --o-alpha-diversity /data/diversity/core_metrics_18000/simpson_vector.qza

echo "Creating alpha-diversity group comparisons..."

for METRIC in shannon simpson evenness
do
  docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:/data" \
    --workdir /data \
    "$QIIME_IMAGE" \
    qiime diversity alpha-group-significance \
      --i-alpha-diversity "/data/diversity/core_metrics_18000/${METRIC}_vector.qza" \
      --m-metadata-file /data/diversity/metadata_biological.tsv \
      --o-visualization "/data/diversity/alpha_significance/${METRIC}_group_significance.qzv"
done

echo "Alpha-diversity analysis finished."
