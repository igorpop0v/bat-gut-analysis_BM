#!/usr/bin/env bash

set -e
set -o pipefail

ASSEMBLY_DIR="$1"
METABULI_DIR="$2"
OUT_DIR="$3"

mkdir -p "$OUT_DIR"

TABLE="$OUT_DIR/parvoviridae_candidates.tsv"
ALL_FASTA="$OUT_DIR/parvoviridae_candidates.fasta"

printf "sample_id\tcontig_id\tlength\ttaxid\tscore\te_value\trank\tlineage\n" > "$TABLE"
: > "$ALL_FASTA"

for CLASSIFICATIONS in \
    "$METABULI_DIR"/RNA_S12046Nr*/RNA_S12046Nr*_classifications.tsv
do
    SAMPLE=$(basename "$(dirname "$CLASSIFICATIONS")")
    CONTIGS="$ASSEMBLY_DIR/$SAMPLE/final.contigs.fa"

    IDS="$OUT_DIR/${SAMPLE}_parvoviridae_ids.txt"
    SAMPLE_FASTA="$OUT_DIR/${SAMPLE}_parvoviridae.fasta"

    echo "Extracting Parvoviridae contigs from $SAMPLE..."

    awk -F '\t' \
      '!/^#/ && $8 ~ /f_Parvoviridae/ {print $2}' \
      "$CLASSIFICATIONS" \
      > "$IDS"

    awk -F '\t' -v sample="$SAMPLE" \
      'BEGIN {OFS="\t"}
       !/^#/ && $8 ~ /f_Parvoviridae/ {
           print sample, $2, $4, $3, $5, $6, $7, $8
       }' \
      "$CLASSIFICATIONS" \
      >> "$TABLE"

    awk -v sample="$SAMPLE" '
        NR == FNR {
            wanted[$1] = 1
            next
        }

        /^>/ {
            header = substr($0, 2)
            split(header, fields, /[ \t]/)
            keep = fields[1] in wanted

            if (keep) {
                print ">" sample "|" header
            }

            next
        }

        keep {
            print
        }
    ' \
      "$IDS" \
      "$CONTIGS" \
      > "$SAMPLE_FASTA"

    cat "$SAMPLE_FASTA" >> "$ALL_FASTA"
done

echo "Candidate extraction completed."
echo "Table: $TABLE"
echo "FASTA: $ALL_FASTA"
