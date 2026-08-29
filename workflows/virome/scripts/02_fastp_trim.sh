#!/usr/bin/env bash

set -e

RAW_DIR="$1"
OUT_DIR="$2"
THREADS="${3:-8}"

mkdir -p "$OUT_DIR/trimmed"
mkdir -p "$OUT_DIR/reports"
mkdir -p "$OUT_DIR/multiqc"

for R1 in "$RAW_DIR"/*.1.fastq.gz
do
  SAMPLE=$(basename "$R1" .1.fastq.gz)
  R2="$RAW_DIR/$SAMPLE.2.fastq.gz"

  echo "Processing $SAMPLE..."

  docker run --rm \
    --user "$(id -u):$(id -g)" \
    --volume "$RAW_DIR:/reads:ro" \
    --volume "$OUT_DIR:/output" \
    quay.io/biocontainers/fastp:1.3.6--h43da1c4_0 \
    fastp \
      --in1 "/reads/$SAMPLE.1.fastq.gz" \
      --in2 "/reads/$SAMPLE.2.fastq.gz" \
      --out1 "/output/trimmed/$SAMPLE.1.fastq.gz" \
      --out2 "/output/trimmed/$SAMPLE.2.fastq.gz" \
      --detect_adapter_for_pe \
      --qualified_quality_phred 20 \
      --length_required 50 \
      --thread "$THREADS" \
      --html "/output/reports/$SAMPLE.fastp.html" \
      --json "/output/reports/$SAMPLE.fastp.json"
done

echo "Creating fastp MultiQC report..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$OUT_DIR:/data" \
  multiqc/multiqc:v1.35 \
  multiqc /data/reports \
    --outdir /data/multiqc \
    --filename fastp_multiqc_report.html \
    --force

echo "fastp completed."
echo "Clean reads: $OUT_DIR/trimmed"
echo "MultiQC report: $OUT_DIR/multiqc/fastp_multiqc_report.html"
