#!/usr/bin/env bash

set -euo pipefail

CURRENT_FASTA="$1"
PREVIOUS_FASTA="$2"
REFSEQ_FASTA="$3"
SOURCE_MANIFEST="$4"
OUTPUT_DIR="$5"

mkdir -p "$OUTPUT_DIR"

OUTPUT_FASTA="${OUTPUT_DIR}/all_densoviruses_final.fasta"
OUTPUT_MANIFEST="${OUTPUT_DIR}/densovirus_final_manifest.tsv"

# Combine the corrected current-study sequences,
# previous-study sequences and the RefSeq genome.
awk -v refseq_file="$REFSEQ_FASTA" '
FILENAME == refseq_file && /^>/ {
    print ">RefSeq_NC_031450.1"
    next
}

{
    print
}
' \
"$CURRENT_FASTA" \
"$PREVIOUS_FASTA" \
"$REFSEQ_FASTA" \
> "$OUTPUT_FASTA"

# Update the length and notes for the corrected Current_Nr11 sequence.
awk -F $'\t' '
BEGIN {
    OFS = "\t"
}

NR == 1 {
    print
    next
}

$1 == "Current_Nr11" {
    $7 = 4904
    $9 = $9 "; read-supported deletion at original position 1997"
}

{
    print
}
' "$SOURCE_MANIFEST" > "$OUTPUT_MANIFEST"

echo "Final densovirus genome set prepared."
echo "Sequences: $(grep -c "^>" "$OUTPUT_FASTA")"
echo "FASTA: $OUTPUT_FASTA"
echo "Manifest: $OUTPUT_MANIFEST"
