#!/usr/bin/env bash

# Remove the positive control, build a phylogenetic tree,
# and generate alpha-rarefaction curves.

set -e

WORK_DIR="$1"
QIIME_IMAGE="quay.io/qiime2/qiime2:2026.7"
DIVERSITY_DIR="$WORK_DIR/diversity"

mkdir -p "$DIVERSITY_DIR"

echo "Removing the positive control..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime feature-table filter-samples \
    --i-table /data/table_deblur_190.qza \
    --m-metadata-file /data/metadata_16s.tsv \
    --p-where "[group] != 'positive_control'" \
    --o-filtered-table /data/diversity/table_biological.qza

echo "Summarizing the biological-sample table..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime feature-table summarize \
    --i-table /data/diversity/table_biological.qza \
    --m-metadata-file /data/metadata_16s.tsv \
    --o-feature-frequencies /data/diversity/feature_frequencies_biological.qza \
    --o-sample-frequencies /data/diversity/sample_frequencies_biological.qza \
    --o-summary /data/diversity/table_biological_summary.qzv

echo "Building the phylogenetic tree..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences /data/rep_seqs_deblur_190.qza \
    --o-alignment /data/diversity/aligned_rep_seqs.qza \
    --o-masked-alignment /data/diversity/masked_aligned_rep_seqs.qza \
    --o-tree /data/diversity/unrooted_tree.qza \
    --o-rooted-tree /data/diversity/rooted_tree.qza

echo "Generating alpha-rarefaction curves..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --volume "$WORK_DIR:/data" \
  --workdir /data \
  "$QIIME_IMAGE" \
  qiime diversity alpha-rarefaction \
    --i-table /data/diversity/table_biological.qza \
    --i-phylogeny /data/diversity/rooted_tree.qza \
    --p-max-depth 18000 \
    --p-min-depth 1000 \
    --p-steps 10 \
    --p-iterations 10 \
    --p-metrics observed_features \
    --p-metrics shannon \
    --p-metrics simpson \
    --p-metrics pielou_e \
    --p-metrics faith_pd \
    --m-metadata-file /data/metadata_16s.tsv \
    --o-visualization /data/diversity/alpha_rarefaction_18000.qzv

echo "Diversity preparation finished."
