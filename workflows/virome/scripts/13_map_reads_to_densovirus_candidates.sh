#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

READS_DIR="$1"
CANDIDATES_FASTA="$2"
OUT_DIR="$3"
THREADS="${4:-12}"

IMAGE="bat-gut-densovirus-mapping:1.0"

mkdir -p \
  "$OUT_DIR/references" \
  "$OUT_DIR/indexes" \
  "$OUT_DIR/bam" \
  "$OUT_DIR/depth" \
  "$OUT_DIR/coverage" \
  "$OUT_DIR/reports"

SUMMARY_FILE="$OUT_DIR/densovirus_mapping_summary.tsv"

printf \
  "candidate_id\tsample_id\tcandidate_length\tmapped_read_segments\tcovered_bases\tcoverage_pct\tmean_depth\tmean_base_quality\tmean_mapping_quality\n" \
  > "$SUMMARY_FILE"

docker run --rm \
  "$IMAGE" \
  bowtie2 --version \
  > "$OUT_DIR/software_versions.txt"

docker run --rm \
  "$IMAGE" \
  samtools --version \
  >> "$OUT_DIR/software_versions.txt"

grep '^>' "$CANDIDATES_FASTA" \
| sed 's/^>//; s/[[:space:]].*$//' \
| while IFS= read -r CANDIDATE_ID
do
    SAMPLE_ID="${CANDIDATE_ID%%|*}"
    SAFE_ID="${CANDIDATE_ID//|/__}"

    R1="$READS_DIR/${SAMPLE_ID}_fwd.fq.gz"
    R2="$READS_DIR/${SAMPLE_ID}_rev.fq.gz"

    CANDIDATE_REFERENCE="$OUT_DIR/references/${SAFE_ID}.fasta"
    INDEX_DIR="$OUT_DIR/indexes/${SAFE_ID}"

    mkdir -p "$INDEX_DIR"

    awk -v target="$CANDIDATE_ID" '
        /^>/ {
            current_id = $0
            sub(/^>/, "", current_id)
            sub(/[[:space:]].*$/, "", current_id)

            keep = current_id == target
        }

        keep {
            print
        }
    ' "$CANDIDATES_FASTA" \
    > "$CANDIDATE_REFERENCE"

    echo "Mapping reads to $CANDIDATE_ID..."

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      --volume "$READS_DIR:/reads:ro" \
      --volume "$OUT_DIR:/output" \
      "$IMAGE" \
      bash -c "
        set -euo pipefail

        bowtie2-build \
          /output/references/${SAFE_ID}.fasta \
          /output/indexes/${SAFE_ID}/${SAFE_ID} \
          > /output/reports/${SAFE_ID}.bowtie2_build.log \
          2>&1

        bowtie2 \
          --very-sensitive-local \
          --no-unal \
          --threads ${THREADS} \
          -x /output/indexes/${SAFE_ID}/${SAFE_ID} \
          -1 /reads/${SAMPLE_ID}_fwd.fq.gz \
          -2 /reads/${SAMPLE_ID}_rev.fq.gz \
          2> /output/reports/${SAFE_ID}.bowtie2.log \
        | samtools view \
            -b \
            -F 4 \
            - \
        | samtools sort \
            -@ ${THREADS} \
            -m 512M \
            -T /output/bam/${SAFE_ID}.temporary \
            -o /output/bam/${SAFE_ID}.bam \
            -

        samtools index \
          /output/bam/${SAFE_ID}.bam

        samtools coverage \
          /output/bam/${SAFE_ID}.bam \
          > /output/coverage/${SAFE_ID}.coverage.tsv

        samtools depth \
          -aa \
          /output/bam/${SAFE_ID}.bam \
          > /output/depth/${SAFE_ID}.depth.tsv

        samtools flagstat \
          -@ ${THREADS} \
          /output/bam/${SAFE_ID}.bam \
          > /output/reports/${SAFE_ID}.flagstat.txt
      "

    COVERAGE_ROW=$(
      awk '
        NR == 2 {
            print
            exit
        }
      ' "$OUT_DIR/coverage/${SAFE_ID}.coverage.tsv"
    )

    IFS=$'\t' read -r \
      CONTIG_NAME \
      START_POSITION \
      END_POSITION \
      MAPPED_READS \
      COVERED_BASES \
      COVERAGE_PERCENT \
      MEAN_DEPTH \
      MEAN_BASE_QUALITY \
      MEAN_MAPPING_QUALITY \
      <<< "$COVERAGE_ROW"

    printf \
      "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$CANDIDATE_ID" \
      "$SAMPLE_ID" \
      "$END_POSITION" \
      "$MAPPED_READS" \
      "$COVERED_BASES" \
      "$COVERAGE_PERCENT" \
      "$MEAN_DEPTH" \
      "$MEAN_BASE_QUALITY" \
      "$MEAN_MAPPING_QUALITY" \
      >> "$SUMMARY_FILE"
done

echo "Densovirus read mapping completed."
echo "Summary: $SUMMARY_FILE"
