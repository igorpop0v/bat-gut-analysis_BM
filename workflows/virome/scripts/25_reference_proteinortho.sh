#!/usr/bin/env bash

set -euo pipefail

PROTEIN_FASTA="$1"
OUT_DIR="$2"
THREADS="${3:-12}"

IMAGE="quay.io/biocontainers/proteinortho:6.3.6--h2b77389_0"

PROTEOME_DIR="$OUT_DIR/proteomes"

mkdir -p "$PROTEOME_DIR"

echo "Creating one protein FASTA per reference genome..."

awk -v output_dir="$PROTEOME_DIR" '
/^>/ {
    header = substr($0, 2)
    split(header, fields, "__")
    output_file = output_dir "/" fields[1] ".faa"
}
{
    print > output_file
}
' "$PROTEIN_FASTA"

PROTEOME_COUNT=$(
    find "$PROTEOME_DIR" \
    -maxdepth 1 \
    -type f \
    -name "*.faa" \
    | wc -l
)

echo "Reference proteomes: $PROTEOME_COUNT"
echo "Running Proteinortho..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --volume "$OUT_DIR:/work" \
    --workdir /work \
    "$IMAGE" \
    bash -lc \
    "proteinortho6.pl proteomes/*.faa \
    -project=densovirus_reference \
    -cpus=$THREADS"

echo "Proteinortho analysis completed."
echo "Results: $OUT_DIR"
