#!/usr/bin/env bash

# Apply QIIME 2 quality filtering before Deblur denoising.

set -e

WORK_DIR="$1"
THREADS="${2:-8}"

QIIME_IMAGE="quay.io/qiime2/qiime2:2026.7"

INPUT="$WORK_DIR/demux_forward_trimmed.qza"
FILTERED="$WORK_DIR/demux_forward_filtered.qza"
FILTERED_SUMMARY="$WORK_DIR/demux_forward_filtered.qzv"

STATS="$WORK_DIR/quality_filter_stats.qza"
STATS_SUMMARY="$WORK_DIR/quality_filter_stats.qzv"

echo "Filtering reads by quality score..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime quality-filter q-score \
    --i-demux "$INPUT" \
    --p-num-processes "$THREADS" \
    --o-filtered-sequences "$FILTERED" \
    --o-filter-stats "$STATS"

echo "Creating filtered-read summary..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime demux summarize \
    --i-data "$FILTERED" \
    --o-visualization "$FILTERED_SUMMARY"

echo "Creating quality-filter statistics table..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime metadata tabulate \
    --m-input-file "$STATS" \
    --o-visualization "$STATS_SUMMARY"

echo "Done."
echo "Filtered reads: $FILTERED"
echo "Filter summary: $FILTERED_SUMMARY"
echo "Filter stats:   $STATS_SUMMARY"
