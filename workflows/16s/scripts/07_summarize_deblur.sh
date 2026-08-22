#!/usr/bin/env bash

# Create summaries of the Deblur feature table and representative sequences.

set -e

WORK_DIR="$1"

QIIME_IMAGE="quay.io/qiime2/qiime2:2026.7"

TABLE="$WORK_DIR/table_deblur_190.qza"
REP_SEQS="$WORK_DIR/rep_seqs_deblur_190.qza"

TABLE_SUMMARY="$WORK_DIR/table_deblur_190_summary.qzv"
FEATURE_FREQUENCIES="$WORK_DIR/feature_frequencies_deblur_190.qza"
SAMPLE_FREQUENCIES="$WORK_DIR/sample_frequencies_deblur_190.qza"
REP_SEQS_SUMMARY="$WORK_DIR/rep_seqs_deblur_190_summary.qzv"

echo "Creating feature table summary..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:$WORK_DIR" \
    --workdir "$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime feature-table summarize \
    --i-table "$TABLE" \
    --o-feature-frequencies "$FEATURE_FREQUENCIES" \
    --o-sample-frequencies "$SAMPLE_FREQUENCIES" \
    --o-summary "$TABLE_SUMMARY"

echo "Creating representative sequence summary..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:$WORK_DIR" \
    --workdir "$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime feature-table tabulate-seqs \
    --i-data "$REP_SEQS" \
    --o-visualization "$REP_SEQS_SUMMARY"

echo "Done."
echo "Feature table summary: $TABLE_SUMMARY"
echo "Sequence summary:      $REP_SEQS_SUMMARY"
