#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

COMBINED_FASTA="$1"
MANIFEST="$2"
OUT_DIR="$3"
THREADS="${4:-12}"

IMAGE="ncbi/blast:2.17.0"

SCRIPT_DIR=$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd
)

PLOT_SCRIPT="$SCRIPT_DIR/../r/11_plot_densovirus_structure.R"

mkdir -p "$OUT_DIR"

REFERENCE_FASTA="$OUT_DIR/previous_p5_reference.fasta"
BLAST_RESULTS="$OUT_DIR/densovirus_vs_previous_p5_hsps.tsv"

: > "$REFERENCE_FASTA"

awk '
  /^>/ {
      current_id = $0
      sub(/^>/, "", current_id)
      sub(/[[:space:]].*$/, "", current_id)

      keep = current_id == "Previous_P5"
  }

  keep {
      print
  }
' "$COMBINED_FASTA" \
> "$REFERENCE_FASTA"

echo "Creating Previous_P5 BLAST database..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$OUT_DIR:/output" \
  "$IMAGE" \
  makeblastdb \
  -in /output/previous_p5_reference.fasta \
  -dbtype nucl \
  -parse_seqids \
  -out /output/previous_p5_db

COMBINED_DIR=$(dirname "$COMBINED_FASTA")
COMBINED_NAME=$(basename "$COMBINED_FASTA")

echo "Comparing all sequences with Previous_P5..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$COMBINED_DIR:/input:ro" \
  --volume "$OUT_DIR:/output" \
  "$IMAGE" \
  blastn \
  -task blastn \
  -query "/input/$COMBINED_NAME" \
  -db /output/previous_p5_db \
  -evalue 1e-10 \
  -word_size 11 \
  -dust no \
  -max_hsps 100 \
  -num_threads "$THREADS" \
  -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore" \
  -out /output/densovirus_vs_previous_p5_hsps.tsv

echo "Creating structural comparison plots..."

Rscript \
  "$PLOT_SCRIPT" \
  "$BLAST_RESULTS" \
  "$MANIFEST" \
  "$OUT_DIR"

echo "Structural comparison completed."
echo "Results: $OUT_DIR"
