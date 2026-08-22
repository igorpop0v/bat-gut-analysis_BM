#!/usr/bin/env bash

# Denoise primer-trimmed and quality-filtered forward reads with Deblur.

set -e

WORK_DIR="$1"
THREADS="${2:-8}"

QIIME_IMAGE="quay.io/qiime2/qiime2:2026.7"

INPUT="$WORK_DIR/demux_forward_filtered.qza"

TABLE="$WORK_DIR/table_deblur_190.qza"
REP_SEQS="$WORK_DIR/rep_seqs_deblur_190.qza"

STATS="$WORK_DIR/deblur_stats_190.qza"
STATS_SUMMARY="$WORK_DIR/deblur_stats_190.qzv"

echo "Running Deblur..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:$WORK_DIR" \
    --workdir "$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime deblur denoise-16S \
    --i-demultiplexed-seqs "$INPUT" \
    --p-trim-length 190 \
    --p-sample-stats \
    --p-jobs-to-start "$THREADS" \
    --o-table "$TABLE" \
    --o-representative-sequences "$REP_SEQS" \
    --o-stats "$STATS"

echo "Creating Deblur statistics visualization..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime deblur visualize-stats \
    --i-deblur-stats "$STATS" \
    --o-visualization "$STATS_SUMMARY"

echo "Done."
echo "Feature table:            $TABLE"
echo "Representative sequences: $REP_SEQS"
echo "Deblur statistics:        $STATS_SUMMARY"
