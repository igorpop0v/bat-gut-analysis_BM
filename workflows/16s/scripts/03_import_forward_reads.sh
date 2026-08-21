#!/usr/bin/env bash

# Import forward reads into QIIME 2 and summarize them.

set -e

RAW_DIR="$1"
WORK_DIR="$2"

QIIME_IMAGE="quay.io/qiime2/qiime2:2026.7"

MANIFEST="$WORK_DIR/manifest_forward.tsv"
DEMUX="$WORK_DIR/demux_forward.qza"
SUMMARY="$WORK_DIR/demux_forward.qzv"

echo "Importing forward reads..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$RAW_DIR:$RAW_DIR:ro" \
    --volume "$WORK_DIR:$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime tools import \
    --type "SampleData[SequencesWithQuality]" \
    --input-path "$MANIFEST" \
    --input-format SingleEndFastqManifestPhred33V2 \
    --output-path "$DEMUX"

echo "Creating QIIME 2 quality summary..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime demux summarize \
    --i-data "$DEMUX" \
    --o-visualization "$SUMMARY"

echo "Done."
echo "Artifact: $DEMUX"
echo "Summary:  $SUMMARY"
