#!/usr/bin/env bash

# Remove the 341F primers from forward 16S reads.

set -e

WORK_DIR="$1"
THREADS="${2:-8}"

QIIME_IMAGE="quay.io/qiime2/qiime2:2026.7"

INPUT="$WORK_DIR/demux_forward.qza"
TRIMMED="$WORK_DIR/demux_forward_trimmed.qza"
SUMMARY="$WORK_DIR/demux_forward_trimmed.qzv"

STATS="$WORK_DIR/cutadapt_forward_stats.qza"
STATS_SUMMARY="$WORK_DIR/cutadapt_forward_stats.qzv"

echo "Removing 341F primers..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime cutadapt trim-single \
    --i-demultiplexed-sequences "$INPUT" \
    --p-front CCTACGGGDGGCWGCAG \
    --p-front CCTAYGGGGYGCWGCAG \
    --p-match-adapter-wildcards \
    --p-discard-untrimmed \
    --p-cores "$THREADS" \
    --o-trimmed-sequences "$TRIMMED" \
    --o-stats "$STATS"

echo "Creating summary of primer-trimmed reads..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime demux summarize \
    --i-data "$TRIMMED" \
    --o-visualization "$SUMMARY"

echo "Creating Cutadapt statistics table..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime metadata tabulate \
    --m-input-file "$STATS" \
    --o-visualization "$STATS_SUMMARY"

echo "Done."
echo "Trimmed reads:   $TRIMMED"
echo "Quality summary: $SUMMARY"
echo "Cutadapt stats:  $STATS_SUMMARY"
