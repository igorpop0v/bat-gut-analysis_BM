#!/usr/bin/env python3

"""Place current-study densovirus contigs in RefSeq coordinate order."""

import csv
import sys
from collections import defaultdict
from pathlib import Path


if len(sys.argv) != 6:
    print(
        "Usage: 38_relinearize_densovirus_genomes.py "
        "<input.fasta> <manifest.tsv> <blast_hsps.tsv> "
        "<output.fasta> <output_manifest.tsv>"
    )
    sys.exit(1)


input_fasta = Path(sys.argv[1])
manifest_file = Path(sys.argv[2])
blast_file = Path(sys.argv[3])
output_fasta = Path(sys.argv[4])
output_manifest = Path(sys.argv[5])


def read_fasta(path):
    records = {}
    sequence_id = None
    sequence_parts = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if sequence_id is not None:
                    records[sequence_id] = "".join(sequence_parts)

                sequence_id = line[1:].split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line)

    if sequence_id is not None:
        records[sequence_id] = "".join(sequence_parts)

    return records


sequences = read_fasta(input_fasta)

with manifest_file.open() as handle:
    manifest_rows = list(
        csv.DictReader(handle, delimiter="\t")
    )

manifest_by_id = {
    row["analysis_id"]: row
    for row in manifest_rows
}


# BLAST columns:
# qseqid sseqid pident length qstart qend qlen
# sstart send slen evalue bitscore
forward_hsps = defaultdict(list)

with blast_file.open() as handle:
    reader = csv.reader(handle, delimiter="\t")

    for row in reader:
        if len(row) < 12:
            continue

        query_id = row[0]
        query_start = int(row[4])
        query_end = int(row[5])
        reference_start = int(row[7])
        reference_end = int(row[8])
        bitscore = float(row[11])

        if query_start < query_end and reference_start < reference_end:
            forward_hsps[query_id].append(
                {
                    "query_start": query_start,
                    "reference_start": reference_start,
                    "bitscore": bitscore
                }
            )


output_rows = []

output_fasta.parent.mkdir(parents=True, exist_ok=True)
output_manifest.parent.mkdir(parents=True, exist_ok=True)


with output_fasta.open("w") as fasta_handle:
    for sequence_id, sequence in sequences.items():
        metadata = manifest_by_id[sequence_id]

        if metadata["dataset"] == "current_study":
            candidates = [
                hsp
                for hsp in forward_hsps[sequence_id]
                if hsp["reference_start"] <= 20
            ]

            if not candidates:
                raise ValueError(
                    f"No forward HSP near the RefSeq origin "
                    f"for {sequence_id}"
                )

            best_hsp = max(
                candidates,
                key=lambda item: item["bitscore"]
            )

            cut_position = best_hsp["query_start"]
            cut_index = cut_position - 1

            new_sequence = (
                sequence[cut_index:]
                + sequence[:cut_index]
            )

            operation = (
                "rotated_to_refseq_coordinate_order"
                if cut_position != 1
                else "unchanged"
            )

            reference_start = best_hsp["reference_start"]

        else:
            new_sequence = sequence
            cut_position = 1
            reference_start = "NA"
            operation = "unchanged"

        fasta_handle.write(f">{sequence_id}\n")

        for position in range(0, len(new_sequence), 80):
            fasta_handle.write(
                new_sequence[position:position + 80] + "\n"
            )

        output_row = dict(metadata)

        output_row.update(
            {
                "annotation_coordinate_operation": operation,
                "original_base_at_new_position_1": cut_position,
                "refseq_position_at_new_position_1": reference_start
            }
        )

        output_rows.append(output_row)


output_columns = list(manifest_rows[0].keys()) + [
    "annotation_coordinate_operation",
    "original_base_at_new_position_1",
    "refseq_position_at_new_position_1"
]


with output_manifest.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=output_columns,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(output_rows)


rotated_count = sum(
    row["annotation_coordinate_operation"]
    == "rotated_to_refseq_coordinate_order"
    for row in output_rows
)


print("Densovirus coordinate standardization completed.")
print(f"Sequences: {len(output_rows)}")
print(f"Current-study sequences rotated: {rotated_count}")
print(f"FASTA: {output_fasta}")
print(f"Manifest: {output_manifest}")
