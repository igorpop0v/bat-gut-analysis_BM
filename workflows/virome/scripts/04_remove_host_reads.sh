#!/usr/bin/env bash
set -e

READS_DIR="$1"
REFERENCE_DIR="$2"
OUT_DIR="$3"
THREADS="${4:-12}"

mkdir -p "$OUT_DIR/reads"
mkdir -p "$OUT_DIR/reports"

for R1 in "$READS_DIR"/*.1.fastq.gz
do
    SAMPLE=$(basename "$R1" .1.fastq.gz)
    R2="$READS_DIR/$SAMPLE.2.fastq.gz"

    echo "Removing host reads from $SAMPLE..."

    docker run --rm \
        --user "$(id -u):$(id -g)" \
        --volume "$READS_DIR:/reads:ro" \
        --volume "$REFERENCE_DIR:/reference:ro" \
        --volume "$OUT_DIR:/output" \
        quay.io/biocontainers/bowtie2:2.5.5--ha27dd3b_0 \
        bowtie2 \
            --very-sensitive \
            --threads "$THREADS" \
            -x /reference/bowtie2_index/nyctalus_host \
            -1 "/reads/$SAMPLE.1.fastq.gz" \
            -2 "/reads/$SAMPLE.2.fastq.gz" \
            --un-conc-gz "/output/reads/$SAMPLE.%.fastq.gz" \
            -S /dev/null \
            2> "$OUT_DIR/reports/$SAMPLE.bowtie2.log"
done

echo "Host-read removal completed."
echo "Host-depleted reads: $OUT_DIR/reads"
echo "Bowtie2 reports: $OUT_DIR/reports"
