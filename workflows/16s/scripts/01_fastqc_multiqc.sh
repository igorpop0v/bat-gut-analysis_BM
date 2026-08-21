#!/usr/bin/env bash

# Run FastQC for all raw 16S FASTQ files and combine reports with MultiQC.

set -e

RAW_DIR="$1"
QC_DIR="$2"
THREADS="${3:-8}"

FASTQC_IMAGE="quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
MULTIQC_IMAGE="multiqc/multiqc:v1.35"

mkdir -p "$QC_DIR/fastqc"
mkdir -p "$QC_DIR/multiqc"

echo "Running FastQC..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$RAW_DIR:/data:ro" \
    --volume "$QC_DIR/fastqc:/output" \
    "$FASTQC_IMAGE" \
    bash -c \
    "find /data -type f -name '*.fastq.gz' \
    -exec fastqc --threads $THREADS --outdir /output {} +"

echo "Running MultiQC..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$QC_DIR/fastqc:/input:ro" \
    --volume "$QC_DIR/multiqc:/output" \
    "$MULTIQC_IMAGE" \
    multiqc /input \
    --outdir /output \
    --filename multiqc_report.html \
    --force

echo "Done."
echo "Report: $QC_DIR/multiqc/multiqc_report.html"
