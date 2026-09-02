#!/usr/bin/env bash

set -e
set -o pipefail

ASSEMBLY_DIR="$1"
QC_DIR="$2"
THREADS="${3:-12}"

QUAST_IMAGE="staphb/quast:5.3.0"
MULTIQC_IMAGE="multiqc/multiqc:v1.35"

mkdir -p "$QC_DIR"

for CONTIGS in "$ASSEMBLY_DIR"/RNA_S12046Nr*/final.contigs.fa
do
    SAMPLE=$(basename "$(dirname "$CONTIGS")")

    echo "Checking $SAMPLE..."

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      --volume "$ASSEMBLY_DIR:/assembly:ro" \
      --volume "$QC_DIR:/qc" \
      "$QUAST_IMAGE" \
      quast.py \
      "/assembly/$SAMPLE/final.contigs.fa" \
      --labels "$SAMPLE" \
      --output-dir "/qc/$SAMPLE" \
      --threads "$THREADS" \
      --min-contig 500 \
      --no-icarus
done

echo "Creating MultiQC report..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$QC_DIR:/qc" \
  "$MULTIQC_IMAGE" \
  multiqc /qc \
  --outdir /qc/multiqc \
  --force

echo "Assembly QC completed."
echo "Results: $QC_DIR"
