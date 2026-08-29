#!/usr/bin/env bash

set -e

RAW_DIR="$1"
QC_DIR="$2"
THREADS="${3:-8}"

mkdir -p "$QC_DIR/fastqc"
mkdir -p "$QC_DIR/multiqc"

echo "Running FastQC..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$RAW_DIR:/reads:ro" \
  --volume "$QC_DIR:/qc" \
  quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0 \
  sh -c "fastqc \
    --threads $THREADS \
    --outdir /qc/fastqc \
    /reads/*.fastq.gz"

echo "Running MultiQC..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$QC_DIR:/qc" \
  multiqc/multiqc:v1.35 \
  multiqc /qc/fastqc \
    --outdir /qc/multiqc \
    --filename multiqc_report.html \
    --force

echo "QC completed."
echo "MultiQC report: $QC_DIR/multiqc/multiqc_report.html"
