#!/usr/bin/env bash

set -e
set -o pipefail

READS_DIR="$1"
ASSEMBLY_DIR="$2"
THREADS="${3:-12}"

IMAGE="quay.io/biocontainers/megahit:1.2.9--h5ca1c30_6"

mkdir -p "$ASSEMBLY_DIR"
mkdir -p "$ASSEMBLY_DIR/logs"

for R1 in "$READS_DIR"/*_fwd.fq.gz
do
    SAMPLE=$(basename "$R1" _fwd.fq.gz)
    R2="$READS_DIR/${SAMPLE}_rev.fq.gz"
    SAMPLE_OUTPUT="$ASSEMBLY_DIR/$SAMPLE"

    if [[ -f "$SAMPLE_OUTPUT/final.contigs.fa" ]]
    then
        echo "Skipping $SAMPLE: assembly already completed."
        continue
    fi

    if [[ -d "$SAMPLE_OUTPUT" ]]
    then
        echo "Incomplete output directory already exists:"
        echo "$SAMPLE_OUTPUT"
        echo "Move or delete this directory before restarting."
        exit 1
    fi

    echo "Assembling $SAMPLE..."

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      --volume "$READS_DIR:/reads:ro" \
      --volume "$ASSEMBLY_DIR:/assembly" \
      "$IMAGE" \
      megahit \
      --presets meta-sensitive \
      --min-contig-len 500 \
      --memory 0.8 \
      --num-cpu-threads "$THREADS" \
      -1 "/reads/${SAMPLE}_fwd.fq.gz" \
      -2 "/reads/${SAMPLE}_rev.fq.gz" \
      -o "/assembly/$SAMPLE" \
      2>&1 | tee "$ASSEMBLY_DIR/logs/${SAMPLE}.megahit.log"

done

echo "Per-sample assemblies completed."
echo "Results: $ASSEMBLY_DIR"
