#!/usr/bin/env bash

# Create a QIIME 2 manifest containing forward reads only.

set -e

RAW_DIR="$1"
MANIFEST="$2"

mkdir -p "$(dirname "$MANIFEST")"

# Write the manifest header.
printf "sample-id\tabsolute-filepath\n" > "$MANIFEST"

# Add one R1 file per sample.
find "$RAW_DIR" \
    -type f \
    -name "*_L001_R1_001.fastq.gz" \
    | sort \
    | while read -r fastq
do
    sample_id="$(basename "$fastq" "_L001_R1_001.fastq.gz")"
    fastq_path="$(realpath "$fastq")"

    printf "%s\t%s\n" "$sample_id" "$fastq_path"
done >> "$MANIFEST"

echo "Manifest created: $MANIFEST"
