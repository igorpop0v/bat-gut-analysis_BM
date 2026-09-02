#!/usr/bin/env python3

import csv
import sys
from pathlib import Path


if len(sys.argv) != 5:
    print(
        "Usage: 16_standardize_densovirus_orientation.py "
        "<input.fasta> <orientation_summary.tsv> "
        "<output.fasta> <output_manifest.tsv>"
    )
    sys.exit(1)


input_fasta = Path(sys.argv[1])
orientation_file = Path(sys.argv[2])
output_fasta = Path(sys.argv[3])
output_manifest = Path(sys.argv[4])


complement_table = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMLVHBNtgcayrkmlvhbn"
)


def reverse_complement(sequence):
    return sequence.translate(complement_table)[::-1]


def read_fasta(path):
    records = {}
    sequence_id = None
    sequence_parts = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()

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


orientations = {}

with orientation_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        orientations[row["query_id"]] = row["dominant_orientation"]

# Previous_P5 was used as the reference and is absent from the BLAST summary.
orientations["Previous_P5"] = "Forward"

sequences = read_fasta(input_fasta)

missing_orientations = sorted(set(sequences) - set(orientations))

if missing_orientations:
    raise ValueError(
        "No orientation information for: "
        + ", ".join(missing_orientations)
    )


output_fasta.parent.mkdir(parents=True, exist_ok=True)

with output_fasta.open("w") as fasta_handle, \
        output_manifest.open("w", newline="") as manifest_handle:

    fieldnames = [
        "sequence_id",
        "original_orientation",
        "operation",
        "length_nt"
    ]

    writer = csv.DictWriter(
        manifest_handle,
        fieldnames=fieldnames,
        delimiter="\t"
    )

    writer.writeheader()

    for sequence_id, sequence in sequences.items():
        orientation = orientations[sequence_id]

        if orientation == "Reverse":
            standardized_sequence = reverse_complement(sequence)
            operation = "reverse_complemented"
        else:
            standardized_sequence = sequence
            operation = "unchanged"

        fasta_handle.write(f">{sequence_id}\n")

        for position in range(0, len(standardized_sequence), 80):
            fasta_handle.write(
                standardized_sequence[position:position + 80] + "\n"
            )

        writer.writerow(
            {
                "sequence_id": sequence_id,
                "original_orientation": orientation,
                "operation": operation,
                "length_nt": len(standardized_sequence)
            }
        )


print(f"Saved standardized FASTA: {output_fasta}")
print(f"Saved orientation manifest: {output_manifest}")
print(f"Sequences processed: {len(sequences)}")
