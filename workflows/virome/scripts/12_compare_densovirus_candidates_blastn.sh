#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

REFERENCE_DIR="$1"
QUERY_FASTA="$2"
OUT_DIR="$3"
THREADS="${4:-12}"

IMAGE="ncbi/blast:2.17.0"

SCRIPT_DIR=$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd
)

SUMMARY_SCRIPT="$SCRIPT_DIR/../r/09_summarize_densovirus_blast.R"

mkdir -p "$OUT_DIR"

COMBINED_REFERENCE="$OUT_DIR/previous_densoviruses.fasta"
REFERENCE_MANIFEST="$OUT_DIR/previous_densoviruses_manifest.tsv"
BLAST_RESULTS="$OUT_DIR/densovirus_blastn_all_hits.tsv"

: > "$COMBINED_REFERENCE"

printf \
  "reference_id\tfilename\tlength\toriginal_header\n" \
  > "$REFERENCE_MANIFEST"

for FASTA in "$REFERENCE_DIR"/*.fasta
do
    REFERENCE_ID=$(basename "$FASTA" .fasta)

    ORIGINAL_HEADER=$(
        awk '
            /^>/ {
                sub(/^>/, "")
                print
                exit
            }
        ' "$FASTA"
    )

    LENGTH=$(
        awk '
            !/^>/ {
                gsub(/[[:space:]]/, "")
                length_sum += length($0)
            }

            END {
                print length_sum
            }
        ' "$FASTA"
    )

    printf \
      "%s\t%s\t%s\t%s\n" \
      "$REFERENCE_ID" \
      "$(basename "$FASTA")" \
      "$LENGTH" \
      "$ORIGINAL_HEADER" \
      >> "$REFERENCE_MANIFEST"

    awk -v reference_id="$REFERENCE_ID" '
        /^>/ {
            print ">" reference_id
            next
        }

        {
            print toupper($0)
        }
    ' "$FASTA" >> "$COMBINED_REFERENCE"
done

echo "Creating BLAST database..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$OUT_DIR:/output" \
  "$IMAGE" \
  makeblastdb \
  -in /output/previous_densoviruses.fasta \
  -dbtype nucl \
  -parse_seqids \
  -out /output/previous_densoviruses_db

QUERY_DIR=$(dirname "$QUERY_FASTA")
QUERY_NAME=$(basename "$QUERY_FASTA")

echo "Comparing candidates with previous densovirus genomes..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$QUERY_DIR:/query:ro" \
  --volume "$OUT_DIR:/output" \
  "$IMAGE" \
  blastn \
  -task blastn \
  -query "/query/$QUERY_NAME" \
  -db /output/previous_densoviruses_db \
  -evalue 1e-20 \
  -max_target_seqs 10 \
  -num_threads "$THREADS" \
  -dust no \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" \
  -out /output/densovirus_blastn_all_hits.tsv

echo "Combining BLAST alignment segments..."

Rscript \
  "$SUMMARY_SCRIPT" \
  "$BLAST_RESULTS" \
  "$OUT_DIR" \
  "$QUERY_FASTA"

echo "Densovirus comparison completed."
echo "Results: $OUT_DIR"
