#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

GENOMES_FASTA="$1"
OUT_DIR="$2"
THREADS="${3:-12}"

IMAGE="ncbi/blast:2.17.0"

SCRIPT_DIR=$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd
)

SUMMARY_SCRIPT="${SCRIPT_DIR}/../r/09_summarize_densovirus_blast.R"

mkdir -p "$OUT_DIR"

FASTA_DIR=$(dirname "$GENOMES_FASTA")
FASTA_NAME=$(basename "$GENOMES_FASTA")

echo "Creating the densovirus BLASTn database..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$FASTA_DIR:/input:ro" \
  --volume "$OUT_DIR:/output" \
  "$IMAGE" \
  makeblastdb \
  -in "/input/$FASTA_NAME" \
  -dbtype nucl \
  -parse_seqids \
  -out /output/densovirus_all_vs_all_db

echo "Running all-versus-all genome comparison..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$FASTA_DIR:/input:ro" \
  --volume "$OUT_DIR:/output" \
  "$IMAGE" \
  blastn \
  -task blastn \
  -query "/input/$FASTA_NAME" \
  -db /output/densovirus_all_vs_all_db \
  -evalue 1e-20 \
  -word_size 11 \
  -dust no \
  -soft_masking false \
  -max_target_seqs 100 \
  -max_hsps 100 \
  -num_threads "$THREADS" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" \
  -out /output/densovirus_all_vs_all_hsps.tsv

echo "Summarizing BLASTn alignment segments..."

Rscript \
  "$SUMMARY_SCRIPT" \
  "$OUT_DIR/densovirus_all_vs_all_hsps.tsv" \
  "$OUT_DIR" \
  "$GENOMES_FASTA"

echo "All-versus-all genome comparison completed."
echo "Results: $OUT_DIR"
