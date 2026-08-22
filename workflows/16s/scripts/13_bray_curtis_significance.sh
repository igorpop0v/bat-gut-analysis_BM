#!/usr/bin/env bash

# Test differences in Bray-Curtis distances among the three study groups.

set -e

WORK_DIR="$1"
QIIME_IMAGE="quay.io/qiime2/qiime2:2026.7"
OUTPUT_DIR="$WORK_DIR/diversity/bray_curtis_significance"

mkdir -p "$OUTPUT_DIR"

echo "Running Bray-Curtis PERMANOVA..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime diversity beta-group-significance \
    --i-distance-matrix /data/diversity/core_metrics_18000/bray_curtis_distance_matrix.qza \
    --m-metadata-file /data/diversity/metadata_biological.tsv \
    --m-metadata-column group \
    --p-method permanova \
    --p-pairwise \
    --p-permutations 9999 \
    --o-visualization /data/diversity/bray_curtis_significance/bray_curtis_permanova.qzv

echo "Running Bray-Curtis PERMDISP..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime diversity beta-group-significance \
    --i-distance-matrix /data/diversity/core_metrics_18000/bray_curtis_distance_matrix.qza \
    --m-metadata-file /data/diversity/metadata_biological.tsv \
    --m-metadata-column group \
    --p-method permdisp \
    --p-pairwise \
    --p-permutations 9999 \
    --o-visualization /data/diversity/bray_curtis_significance/bray_curtis_permdisp.qzv

echo "Bray-Curtis significance tests finished."

