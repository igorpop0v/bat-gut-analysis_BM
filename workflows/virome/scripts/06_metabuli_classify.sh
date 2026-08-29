#!/usr/bin/env bash

set -e
set -o pipefail

READS_DIR="$1"
DB_DIR="$2"
OUT_DIR="$3"
THREADS="${4:-12}"
MAX_RAM="${5:-48}"

IMAGE="ghcr.io/steineggerlab/metabuli:1.2.0"

mkdir -p "$OUT_DIR"

for R1 in "$READS_DIR"/*_fwd.fq.gz
do
    SAMPLE=$(basename "$R1" _fwd.fq.gz)
    R2="$READS_DIR/${SAMPLE}_rev.fq.gz"
    SAMPLE_OUT="$OUT_DIR/$SAMPLE"

    if [[ -f "$SAMPLE_OUT/${SAMPLE}_report.tsv" ]]
    then
        echo "Skipping $SAMPLE: MetaBuli report already exists."
        continue
    fi

    mkdir -p "$SAMPLE_OUT"

    echo "Classifying $SAMPLE..."

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      --volume "$READS_DIR:/reads:ro" \
      --volume "$DB_DIR:/db:ro" \
      --volume "$OUT_DIR:/output" \
      "$IMAGE" \
      classify \
      "/reads/${SAMPLE}_fwd.fq.gz" \
      "/reads/${SAMPLE}_rev.fq.gz" \
      /db \
      "/output/$SAMPLE" \
      "$SAMPLE" \
      --threads "$THREADS" \
      --max-ram "$MAX_RAM" \
      --precise 1 \
      --validate-input 1 \
      2>&1 | tee "$SAMPLE_OUT/metabuli.log"
done

echo "MetaBuli classification completed."
echo "Results: $OUT_DIR"
