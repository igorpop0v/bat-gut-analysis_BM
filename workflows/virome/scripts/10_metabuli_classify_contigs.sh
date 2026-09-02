#!/usr/bin/env bash

set -e
set -o pipefail

ASSEMBLY_DIR="$1"
DB_DIR="$2"
OUT_DIR="$3"
THREADS="${4:-12}"
MAX_RAM="${5:-48}"

IMAGE="ghcr.io/steineggerlab/metabuli:1.2.0"

mkdir -p "$OUT_DIR"

for CONTIGS in "$ASSEMBLY_DIR"/RNA_S12046Nr*/final.contigs.fa
do
    SAMPLE=$(basename "$(dirname "$CONTIGS")")
    SAMPLE_OUT="$OUT_DIR/$SAMPLE"

    if [[ -f "$SAMPLE_OUT/${SAMPLE}_report.tsv" ]]
    then
        echo "Skipping $SAMPLE: report already exists."
        continue
    fi

    mkdir -p "$SAMPLE_OUT"

    echo "Classifying contigs from $SAMPLE..."

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      --volume "$ASSEMBLY_DIR:/assembly:ro" \
      --volume "$DB_DIR:/db:ro" \
      --volume "$OUT_DIR:/output" \
      "$IMAGE" \
      classify \
      --seq-mode 3 \
      --lineage 1 \
      --threads "$THREADS" \
      --max-ram "$MAX_RAM" \
      --validate-input 1 \
      "/assembly/$SAMPLE/final.contigs.fa" \
      /db \
      "/output/$SAMPLE" \
      "$SAMPLE" \
      2>&1 | tee "$SAMPLE_OUT/metabuli.log"
done

echo "Contig classification completed."
echo "Results: $OUT_DIR"
