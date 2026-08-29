#!/usr/bin/env bash

set -e
set -o pipefail

READS_DIR="$1"
SORTMERNA_DIR="$2"
THREADS="${3:-12}"

IMAGE="bat-gut-sortmerna:7.0.0"

mkdir -p "$SORTMERNA_DIR/non_rRNA"
mkdir -p "$SORTMERNA_DIR/work"
mkdir -p "$SORTMERNA_DIR/logs"

for R1 in "$READS_DIR"/*.1.fastq.gz
do
    SAMPLE=$(basename "$R1" .1.fastq.gz)
    R2="$READS_DIR/$SAMPLE.2.fastq.gz"

    OUT_FWD="$SORTMERNA_DIR/non_rRNA/${SAMPLE}_fwd.fq.gz"
    OUT_REV="$SORTMERNA_DIR/non_rRNA/${SAMPLE}_rev.fq.gz"

    if [[ -f "$OUT_FWD" && -f "$OUT_REV" ]]
    then
        echo "Skipping $SAMPLE: output files already exist."
        continue
    fi

    echo "Removing rRNA from $SAMPLE..."

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      --volume "$READS_DIR:/reads:ro" \
      --volume "$SORTMERNA_DIR:/sortmerna" \
      "$IMAGE" \
      --ref /sortmerna/database/smr_v4.3_default_db.fasta \
      --reads "/reads/$SAMPLE.1.fastq.gz" \
      --reads "/reads/$SAMPLE.2.fastq.gz" \
      --workdir "/sortmerna/work/$SAMPLE" \
      --idx-dir /sortmerna/index \
      --other "/sortmerna/non_rRNA/$SAMPLE" \
      --fastx \
      --paired_in \
      --out2 \
      --threads "$THREADS" \
      2>&1 | tee "$SORTMERNA_DIR/logs/$SAMPLE.sortmerna.log"
done

echo "rRNA removal completed."
echo "Non-rRNA reads: $SORTMERNA_DIR/non_rRNA"
