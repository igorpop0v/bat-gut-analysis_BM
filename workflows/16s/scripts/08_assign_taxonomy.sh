#!/usr/bin/env bash

# Assign GTDB taxonomy to Deblur representative sequences.

set -e

WORK_DIR="$1"
CLASSIFIER="$2"
THREADS="${3:-2}"

CLASSIFIER_DIR="$(dirname "$CLASSIFIER")"

QIIME_IMAGE="quay.io/qiime2/qiime2:2026.7"

REP_SEQS="$WORK_DIR/rep_seqs_deblur_190.qza"

TAXONOMY="$WORK_DIR/taxonomy_gtdb_r232.qza"
TAXONOMY_SUMMARY="$WORK_DIR/taxonomy_gtdb_r232.qzv"

echo "Assigning taxonomy with GTDB r232..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:$WORK_DIR" \
    --volume "$CLASSIFIER_DIR:$CLASSIFIER_DIR:ro" \
    --workdir "$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads "$REP_SEQS" \
    --p-n-jobs "$THREADS" \
    --p-read-orientation auto \
    --o-classification "$TAXONOMY"

echo "Creating taxonomy visualization..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$WORK_DIR:$WORK_DIR" \
    --workdir "$WORK_DIR" \
    "$QIIME_IMAGE" \
    qiime metadata tabulate \
    --m-input-file "$TAXONOMY" \
    --o-visualization "$TAXONOMY_SUMMARY"

echo "Done."
echo "Taxonomy:         $TAXONOMY"
echo "Taxonomy summary: $TAXONOMY_SUMMARY"
