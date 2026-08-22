#!/usr/bin/env bash

# Calculate core alpha and beta diversity metrics at 18,000 reads per sample.

set -e

WORK_DIR="$1"
QIIME_IMAGE="quay.io/qiime2/qiime2:2026.7"
DIVERSITY_DIR="$WORK_DIR/diversity"

mkdir -p "$DIVERSITY_DIR"

# Create metadata containing only biological samples.
grep -v '^Control-16S' \
  "$WORK_DIR/metadata_16s.tsv" \
  > "$DIVERSITY_DIR/metadata_biological.tsv"

echo "Running core diversity metrics at 18,000 reads per sample..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime diversity core-metrics-phylogenetic \
    --i-table /data/diversity/table_biological.qza \
    --i-phylogeny /data/diversity/rooted_tree.qza \
    --p-sampling-depth 18000 \
    --p-random-seed 42 \
    --m-metadata-file /data/diversity/metadata_biological.tsv \
    --output-dir /data/diversity/core_metrics_18000

echo "Core diversity analysis finished."
