#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

CURRENT_CANDIDATES_FASTA="$1"
PREVIOUS_GENOMES_DIR="$2"
OUT_DIR="$3"

SCRIPT_DIR=$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd
)

CONFIG_FILE="$SCRIPT_DIR/../config/densovirus_near_complete_candidates.tsv"

mkdir -p "$OUT_DIR"

NEW_FASTA="$OUT_DIR/current_near_complete_densoviruses.fasta"
PREVIOUS_FASTA="$OUT_DIR/previous_study_densoviruses.fasta"
COMBINED_FASTA="$OUT_DIR/all_densoviruses_combined.fasta"
MANIFEST="$OUT_DIR/densovirus_sequence_manifest.tsv"

: > "$NEW_FASTA"
: > "$PREVIOUS_FASTA"
: > "$COMBINED_FASTA"

printf \
  "analysis_id\tdataset\toriginal_id\tsample_id\tgroup\tstatus\tlength_nt\tsource_file\tnotes\n" \
  > "$MANIFEST"

echo "Extracting current near-complete candidates..."

while IFS=$'\t' read -r \
  CANDIDATE_ID \
  ANALYSIS_ID \
  SAMPLE_ID \
  GROUP \
  STATUS \
  NOTES
do
    if [[ "$CANDIDATE_ID" == "candidate_id" ]]
    then
        continue
    fi

    SEQUENCE_LENGTH=$(
      awk -v target="$CANDIDATE_ID" '
        /^>/ {
            current_id = $0
            sub(/^>/, "", current_id)
            sub(/[[:space:]].*$/, "", current_id)

            keep = current_id == target
            next
        }

        keep {
            gsub(/[[:space:]]/, "")
            sequence_length += length($0)
        }

        END {
            print sequence_length + 0
        }
      ' "$CURRENT_CANDIDATES_FASTA"
    )

    if [[ "$SEQUENCE_LENGTH" -eq 0 ]]
    then
        echo "Candidate not found: $CANDIDATE_ID" >&2
        exit 1
    fi

    for OUTPUT_FASTA in "$NEW_FASTA" "$COMBINED_FASTA"
    do
        printf ">%s\n" "$ANALYSIS_ID" \
          >> "$OUTPUT_FASTA"

        awk -v target="$CANDIDATE_ID" '
          /^>/ {
              current_id = $0
              sub(/^>/, "", current_id)
              sub(/[[:space:]].*$/, "", current_id)

              keep = current_id == target
              next
          }

          keep {
              print toupper($0)
          }
        ' "$CURRENT_CANDIDATES_FASTA" \
          >> "$OUTPUT_FASTA"
    done

    printf \
      "%s\tcurrent_study\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$ANALYSIS_ID" \
      "$CANDIDATE_ID" \
      "$SAMPLE_ID" \
      "$GROUP" \
      "$STATUS" \
      "$SEQUENCE_LENGTH" \
      "$(basename "$CURRENT_CANDIDATES_FASTA")" \
      "$NOTES" \
      >> "$MANIFEST"

done < "$CONFIG_FILE"

echo "Adding genomes from the previous study..."

for FASTA in "$PREVIOUS_GENOMES_DIR"/*.fasta
do
    FILE_ID=$(basename "$FASTA" .fasta)
    SHORT_ID="${FILE_ID#Densovirus_}"
    ANALYSIS_ID="Previous_${SHORT_ID}"

    ORIGINAL_ID=$(
      awk '
        /^>/ {
            sub(/^>/, "")
            sub(/[[:space:]].*$/, "")
            print
            exit
        }
      ' "$FASTA"
    )

    SEQUENCE_LENGTH=$(
      awk '
        !/^>/ {
            gsub(/[[:space:]]/, "")
            sequence_length += length($0)
        }

        END {
            print sequence_length + 0
        }
      ' "$FASTA"
    )

    for OUTPUT_FASTA in "$PREVIOUS_FASTA" "$COMBINED_FASTA"
    do
        printf ">%s\n" "$ANALYSIS_ID" \
          >> "$OUTPUT_FASTA"

        awk '
          !/^>/ {
              print toupper($0)
          }
        ' "$FASTA" \
          >> "$OUTPUT_FASTA"
    done

    printf \
      "%s\tprevious_study\t%s\tNA\tNA\treference_previous_study\t%s\t%s\tGenome from previous study\n" \
      "$ANALYSIS_ID" \
      "$ORIGINAL_ID" \
      "$SEQUENCE_LENGTH" \
      "$(basename "$FASTA")" \
      >> "$MANIFEST"
done

echo "Comparative densovirus dataset prepared."
echo "Current sequences: $NEW_FASTA"
echo "Previous sequences: $PREVIOUS_FASTA"
echo "Combined sequences: $COMBINED_FASTA"
echo "Manifest: $MANIFEST"
