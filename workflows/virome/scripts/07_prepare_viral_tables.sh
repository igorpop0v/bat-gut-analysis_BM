#!/usr/bin/env bash

set -e

export LC_ALL=C

REPORT_DIR="$1"
SORTMERNA_DIR="$2"
OUT_DIR="$3"

mkdir -p "$OUT_DIR"

QC_FILE="$OUT_DIR/sample_qc.tsv"
FAMILY_FILE="$OUT_DIR/viral_family_counts_long.tsv"

printf "sample_id\thost_depleted_pairs\tmetabuli_input_pairs\trrna_removed_pct\tretained_pairs_pct\tclassified_pairs\tclassified_pct\tviral_pairs\tviral_pct_non_rrna\tviral_pct_host_depleted\n" \
  > "$QC_FILE"

printf "sample_id\ttaxid\tfamily\tcount\trelative_abundance_pct\n" \
  > "$FAMILY_FILE"

for REPORT in "$REPORT_DIR"/*/*_report.tsv
do
    SAMPLE=$(basename "$REPORT" _report.tsv)
    SORTMERNA_LOG="$SORTMERNA_DIR/work/$SAMPLE/out/aligned.log"

    HOST_READS=$(awk -F'= ' \
      '/Total reads =/ {print $2}' \
      "$SORTMERNA_LOG")

    HOST_PAIRS=$((HOST_READS / 2))

    RRNA_PCT=$(awk \
      '/Total reads passing E-value threshold/ {
          value=$NF
          gsub(/[()]/, "", value)
          print value
      }' \
      "$SORTMERNA_LOG")

    read -r TOTAL CLASSIFIED VIRAL <<< "$(awk -F'\t' '
        BEGIN {
            unclassified=0
            classified=0
            viral=0
        }

        $5 == 0 {
            unclassified=$2
        }

        $5 == 1 {
            classified=$2
        }

        $5 == 10239 {
            viral=$2
        }

        END {
            total=unclassified + classified
            print total, classified, viral
        }
    ' "$REPORT")"

    RETAINED_PCT=$(awk -v n="$TOTAL" -v d="$HOST_PAIRS" \
      'BEGIN {printf "%.6f", 100*n/d}')

    CLASSIFIED_PCT=$(awk -v n="$CLASSIFIED" -v d="$TOTAL" \
      'BEGIN {printf "%.6f", 100*n/d}')

    VIRAL_PCT=$(awk -v n="$VIRAL" -v d="$TOTAL" \
      'BEGIN {printf "%.6f", 100*n/d}')

    VIRAL_HOST_PCT=$(awk -v n="$VIRAL" -v d="$HOST_PAIRS" \
      'BEGIN {printf "%.6f", 100*n/d}')

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$SAMPLE" \
      "$HOST_PAIRS" \
      "$TOTAL" \
      "$RRNA_PCT" \
      "$RETAINED_PCT" \
      "$CLASSIFIED" \
      "$CLASSIFIED_PCT" \
      "$VIRAL" \
      "$VIRAL_PCT" \
      "$VIRAL_HOST_PCT" \
      >> "$QC_FILE"

    awk -F'\t' -v sample="$SAMPLE" '
        function indentation(text) {
            match(text, /[^ ]/)
            return RSTART - 1
        }

        NR == 1 {
            next
        }

        $5 == 10239 {
            inside_viruses=1
            virus_indent=indentation($6)
            viral_total=$2
            next
        }

        inside_viruses {
            current_indent=indentation($6)

            if (current_indent <= virus_indent) {
                inside_viruses=0
                next
            }

            if ($4 == "family") {
                family_name=$6
                sub(/^ +/, "", family_name)

                family_sum += $2

                printf "%s\t%s\t%s\t%s\t%.6f\n",
                    sample,
                    $5,
                    family_name,
                    $2,
                    100*$2/viral_total
            }
        }

        END {
            unclassified=viral_total-family_sum

            if (unclassified < 0) {
                unclassified=0
            }

            if (viral_total > 0) {
                printf "%s\t0\tUnclassified viruses\t%s\t%.6f\n",
                    sample,
                    unclassified,
                    100*unclassified/viral_total
            }
        }
    ' "$REPORT" >> "$FAMILY_FILE"
done

echo "Viral tables created."
echo "QC table: $QC_FILE"
echo "Family counts: $FAMILY_FILE"
