#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

QUERY_FASTA="$1"
MANIFEST="$2"
REFSEQ_FASTA="$3"
OUT_DIR="$4"
THREADS="${5:-12}"

IMAGE="ncbi/blast:2.17.0"

SCRIPT_DIR=$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd
)

SUMMARY_SCRIPT="$SCRIPT_DIR/../r/12_summarize_refseq_comparison.R"

mkdir -p "$OUT_DIR"

QUERY_DIR=$(dirname "$QUERY_FASTA")
QUERY_NAME=$(basename "$QUERY_FASTA")

REFSEQ_DIR=$(dirname "$REFSEQ_FASTA")
REFSEQ_NAME=$(basename "$REFSEQ_FASTA")

BLAST_RESULTS="$OUT_DIR/densovirus_vs_NC_031450.1_hsps.tsv"

echo "Creating RefSeq BLAST database..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$REFSEQ_DIR:/reference:ro" \
  --volume "$OUT_DIR:/output" \
  "$IMAGE" \
  makeblastdb \
  -in "/reference/$REFSEQ_NAME" \
  -dbtype nucl \
  -parse_seqids \
  -out /output/NC_031450.1_db

echo "Comparing 17 densovirus sequences with NC_031450.1..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$QUERY_DIR:/query:ro" \
  --volume "$OUT_DIR:/output" \
  "$IMAGE" \
  blastn \
  -task blastn \
  -query "/query/$QUERY_NAME" \
  -db /output/NC_031450.1_db \
  -evalue 1e-20 \
  -word_size 11 \
  -dust no \
  -max_hsps 100 \
  -num_threads "$THREADS" \
  -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore" \
  -out /output/densovirus_vs_NC_031450.1_hsps.tsv

echo "Summarizing the comparison..."

Rscript \
  "$SUMMARY_SCRIPT" \
  "$BLAST_RESULTS" \
  "$MANIFEST" \
  "$OUT_DIR"

echo "Comparison completed."
echo "Results: $OUT_DIR"
