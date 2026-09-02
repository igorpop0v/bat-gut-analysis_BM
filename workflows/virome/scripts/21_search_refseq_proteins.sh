#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

PROTEIN_FASTA="$1"
PROTEIN_MANIFEST="$2"
GENOME_FASTA="$3"
GENOME_MANIFEST="$4"
OUT_DIR="$5"
THREADS="${6:-12}"

IMAGE="ncbi/blast:2.17.0"

SCRIPT_DIR=$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd
)

SUMMARY_SCRIPT="$SCRIPT_DIR/20_summarize_refseq_protein_hits.py"

mkdir -p "$OUT_DIR"

PROTEIN_DIR=$(dirname "$PROTEIN_FASTA")
PROTEIN_NAME=$(basename "$PROTEIN_FASTA")

GENOME_DIR=$(dirname "$GENOME_FASTA")
GENOME_NAME=$(basename "$GENOME_FASTA")

RAW_RESULTS="$OUT_DIR/refseq_proteins_tblastn_all_hits.tsv"
LONG_RESULTS="$OUT_DIR/refseq_protein_hits_by_genome.tsv"
SUMMARY_RESULTS="$OUT_DIR/refseq_marker_completeness_summary.tsv"

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
  -out /output/densovirus_genomes_db

GENOME_COUNT=$(grep -c '^>' "$GENOME_FASTA")
echo "Searching RefSeq proteins in ${GENOME_COUNT} genomes..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$PROTEIN_DIR:/proteins:ro" \
  --volume "$OUT_DIR:/output" \
  "$IMAGE" \
  tblastn \
  -query "/proteins/$PROTEIN_NAME" \
  -db /output/densovirus_genomes_db \
  -evalue 1e-5 \
  -seg no \
  -max_hsps 20 \
  -num_threads "$THREADS" \
  -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore sframe" \
  -out /output/refseq_proteins_tblastn_all_hits.tsv

echo "Summarizing protein hits..."

python3 \
  "$SUMMARY_SCRIPT" \
  "$RAW_RESULTS" \
  "$PROTEIN_MANIFEST" \
  "$GENOME_MANIFEST" \
  "$LONG_RESULTS" \
  "$SUMMARY_RESULTS"

echo "Protein search completed."
echo "Results: $OUT_DIR"
