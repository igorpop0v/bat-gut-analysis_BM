#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

QUERY_PROTEINS="$1"
REFERENCE_PROTEINS="$2"
OUT_DIR="$3"
THREADS="${4:-12}"

IMAGE="ncbi/blast:2.17.0"

mkdir -p "$OUT_DIR"

QUERY_DIR=$(dirname "$QUERY_PROTEINS")
QUERY_NAME=$(basename "$QUERY_PROTEINS")
REFERENCE_DIR=$(dirname "$REFERENCE_PROTEINS")
REFERENCE_NAME=$(basename "$REFERENCE_PROTEINS")

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$REFERENCE_DIR:/reference:ro" \
  --volume "$OUT_DIR:/output" \
  "$IMAGE" \
  makeblastdb \
  -in "/reference/$REFERENCE_NAME" \
  -dbtype prot \
  -parse_seqids \
  -out /output/refseq_densovirus_proteins

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "$QUERY_DIR:/query:ro" \
  --volume "$OUT_DIR:/output" \
  "$IMAGE" \
  blastp \
  -query "/query/$QUERY_NAME" \
  -db /output/refseq_densovirus_proteins \
  -evalue 1e-10 \
  -seg no \
  -max_target_seqs 5 \
  -num_threads "$THREADS" \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
  -out /output/provisional_proteins_blastp_all_hits.tsv

awk -F '\t' '
BEGIN {
  OFS = "\t"
  print "query_id", "best_reference_id", "expected_orf", "reference_orf", \
        "identity_pct", "query_coverage_pct", "reference_coverage_pct", \
        "e_value", "bitscore", "validation_status"
}

!seen[$1]++ {
  split($1, query_parts, "_")
  expected_orf = query_parts[length(query_parts)]

  split($2, reference_parts, "|")
  split(reference_parts[2], alias_parts, "_")
  reference_orf = "NA"
  for (i = 1; i <= length(alias_parts); i++) {
    if (alias_parts[i] ~ /^ORF[1-5]$/) {
      reference_orf = alias_parts[i]
    }
  }

  query_coverage = 100 * $4 / $5
  reference_coverage = 100 * $4 / $6

  status = "REVIEW"
  if (expected_orf == reference_orf && $3 >= 90 && query_coverage >= 95 && reference_coverage >= 95) {
    status = "PASS"
  }

  print $1, $2, expected_orf, reference_orf, $3, \
        sprintf("%.3f", query_coverage), \
        sprintf("%.3f", reference_coverage), $7, $8, status
}
' \
"$OUT_DIR/provisional_proteins_blastp_all_hits.tsv" \
> "$OUT_DIR/provisional_proteins_blastp_best_hits.tsv"

QUERY_COUNT=$(grep -c '^>' "$QUERY_PROTEINS")
PASS_COUNT=$(awk -F '\t' '$10 == "PASS" {count++} END {print count + 0}' \
  "$OUT_DIR/provisional_proteins_blastp_best_hits.tsv")

echo "Protein validation completed."
echo "Proteins tested: $QUERY_COUNT"
echo "Proteins passing validation: $PASS_COUNT"
echo "Results: $OUT_DIR/provisional_proteins_blastp_best_hits.tsv"
