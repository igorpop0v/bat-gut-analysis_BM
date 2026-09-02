#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

PREVIOUS_LIST="$1"
OUT_DIR="$2"

IMAGE="ncbi/edirect:latest"

mkdir -p "$OUT_DIR"

PREVIOUS_SORTED="$OUT_DIR/refseq_accessions_2024.txt"
CURRENT_LIST="$OUT_DIR/refseq_accessions_current.txt"
REFSEQ_UNION="$OUT_DIR/refseq_accessions_union.txt"
FINAL_PANEL="$OUT_DIR/reference_panel_accessions.txt"

NEW_ACCESSIONS="$OUT_DIR/refseq_accessions_added_to_2024_panel.txt"
OLD_OUTSIDE_CURRENT="$OUT_DIR/refseq_accessions_outside_current_query.txt"

GENOME_FASTA="$OUT_DIR/reference_panel_genomes.fasta"
GENBANK_FILE="$OUT_DIR/reference_panel_records.gb"

QUERY='txid40120[Organism:exp] AND "complete genome"[Title] AND viruses[filter] AND refseq[filter] AND 3000:6500[SLEN]'

echo "Preparing the previous RefSeq accession list..."

sort -u "$PREVIOUS_LIST" > "$PREVIOUS_SORTED"

echo "Searching the current RefSeq Densovirinae collection..."

docker run --rm \
  --env QUERY="$QUERY" \
  "$IMAGE" \
  /bin/sh -c \
  'esearch -db nucleotide -query "$QUERY" | efetch -format acc' \
  | sed '/^[[:space:]]*$/d' \
  | sort -u \
  > "$CURRENT_LIST"

echo "Combining previous and current RefSeq accessions..."

cat \
  "$PREVIOUS_SORTED" \
  "$CURRENT_LIST" \
  | sort -u \
  > "$REFSEQ_UNION"

comm -13 \
  "$PREVIOUS_SORTED" \
  "$CURRENT_LIST" \
  > "$NEW_ACCESSIONS"

comm -23 \
  "$PREVIOUS_SORTED" \
  "$CURRENT_LIST" \
  > "$OLD_OUTSIDE_CURRENT"

echo "Adding the Tenebrio molitor densovirus genome..."

{
  cat "$REFSEQ_UNION"
  echo "MW628494.1"
} | sort -u > "$FINAL_PANEL"

echo "Downloading nucleotide sequences..."

docker run --rm -i \
  "$IMAGE" \
  /bin/sh -c \
  'epost -db nucleotide | efetch -format fasta' \
  < "$FINAL_PANEL" \
  > "$GENOME_FASTA"

echo "Downloading GenBank records..."

docker run --rm -i \
  "$IMAGE" \
  /bin/sh -c \
  'epost -db nucleotide | efetch -format gbwithparts' \
  < "$FINAL_PANEL" \
  > "$GENBANK_FILE"

ACCESSION_COUNT=$(wc -l < "$FINAL_PANEL")
FASTA_COUNT=$(grep -c '^>' "$GENOME_FASTA")

echo
echo "Previous RefSeq accessions: $(wc -l < "$PREVIOUS_SORTED")"
echo "Current Densovirinae RefSeq accessions: $(wc -l < "$CURRENT_LIST")"
echo "Unique RefSeq accessions: $(wc -l < "$REFSEQ_UNION")"
echo "Final external reference panel: $ACCESSION_COUNT"
echo "Downloaded FASTA sequences: $FASTA_COUNT"

if [ "$ACCESSION_COUNT" -ne "$FASTA_COUNT" ]; then
  echo "Warning: accession and FASTA sequence counts differ."
  exit 1
fi

echo
echo "Reference panel preparation completed."
echo "Results: $OUT_DIR"
