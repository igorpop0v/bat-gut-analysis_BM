#!/usr/bin/env bash

set -euo pipefail

REFERENCE_PROTEINS="$1"
GENOME_FASTA="$2"
GENOME_MANIFEST="$3"
OUT_DIR="$4"
THREADS="${5:-12}"

IMAGE="ncbi/blast:2.17.0"

SCRIPT_DIR=$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd
)

PARSER="$SCRIPT_DIR/28_reconstruct_study_marker_proteins.py"

mkdir -p "$OUT_DIR"

REFERENCE_DIR=$(dirname "$REFERENCE_PROTEINS")
REFERENCE_NAME=$(basename "$REFERENCE_PROTEINS")

GENOME_DIR=$(dirname "$GENOME_FASTA")
GENOME_NAME=$(basename "$GENOME_FASTA")

QUERY_FASTA="$OUT_DIR/nsp1_sp1_queries.faa"
RAW_RESULTS="$OUT_DIR/study_markers_tblastn.tsv"

awk '
/^>/ {keep = 0}
/^>YP_009310053.1[| ]/ {keep = 1}
/^>YP_009310055.1[| ]/ {keep = 1}
keep
' "$REFERENCE_PROTEINS" > "$QUERY_FASTA"

QUERY_COUNT=$(grep -c '^>' "$QUERY_FASTA")

if [[ "$QUERY_COUNT" -ne 2 ]]; then
  echo "Expected 2 marker queries, found $QUERY_COUNT" >&2
  exit 1
fi

echo "Creating nucleotide BLAST database..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$GENOME_DIR:/genomes:ro" \
  --volume "$OUT_DIR:/output" \
  "$IMAGE" \
  makeblastdb \
  -in "/genomes/$GENOME_NAME" \
  -dbtype nucl \
  -parse_seqids \
  -out /output/study_genomes_db

echo "Searching for NSP1 and capsid proteins..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$OUT_DIR:/output" \
  "$IMAGE" \
  tblastn \
  -query /output/nsp1_sp1_queries.faa \
  -db /output/study_genomes_db \
  -evalue 1e-10 \
  -seg no \
  -max_hsps 20 \
  -num_threads "$THREADS" \
  -outfmt "6 qseqid sseqid qstart qend qlen sstart send sframe pident evalue bitscore qseq sseq" \
  -out /output/study_markers_tblastn.tsv

echo "Reconstructing marker proteins..."

python3 \
  "$PARSER" \
  "$RAW_RESULTS" \
  "$GENOME_MANIFEST" \
  "$OUT_DIR"

echo "Study marker extraction completed."
echo "Results: $OUT_DIR"
