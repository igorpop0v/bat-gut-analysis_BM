#!/usr/bin/env python3

import csv
import statistics
import sys
from pathlib import Path


input_dir = Path(sys.argv[1])
output_file = Path(sys.argv[2])


def read_fasta(fasta_file):

    sequences = []
    current_sequence = []

    with open(fasta_file) as handle:

        for line in handle:

            line = line.strip()

            if line.startswith(">"):

                if current_sequence:
                    sequences.append("".join(current_sequence))
                    current_sequence = []

            elif line:
                current_sequence.append(line)

    if current_sequence:
        sequences.append("".join(current_sequence))

    return sequences


rows = []

for fasta_file in sorted(input_dir.glob("*.faa")):

    sequences = read_fasta(fasta_file)

    lengths = {len(sequence) for sequence in sequences}

    if len(lengths) != 1:
        raise ValueError(
            f"Sequences have different aligned lengths in {fasta_file}"
        )

    alignment_length = lengths.pop()
    total_positions = len(sequences) * alignment_length

    gap_count = sum(sequence.count("-") for sequence in sequences)

    unknown_count = sum(
        sequence.upper().count("X")
        for sequence in sequences
    )

    ungapped_lengths = [
        len(sequence.replace("-", ""))
        for sequence in sequences
    ]

    column_occupancy = []

    for position in range(alignment_length):

        present = sum(
            sequence[position] not in {"-", "X", "x", "?"}
            for sequence in sequences
        )

        column_occupancy.append(
            present / len(sequences)
        )

    rows.append(
        {
            "alignment": fasta_file.name,
            "sequences": len(sequences),
            "alignment_length": alignment_length,
            "gap_pct": round(
                100 * gap_count / total_positions,
                3
            ),
            "unknown_pct": round(
                100 * unknown_count / total_positions,
                3
            ),
            "minimum_ungapped_length": min(ungapped_lengths),
            "median_ungapped_length": round(
                statistics.median(ungapped_lengths),
                1
            ),
            "maximum_ungapped_length": max(ungapped_lengths),
            "columns_ge_50_pct_occupancy": sum(
                value >= 0.50
                for value in column_occupancy
            ),
            "columns_ge_70_pct_occupancy": sum(
                value >= 0.70
                for value in column_occupancy
            ),
            "columns_ge_90_pct_occupancy": sum(
                value >= 0.90
                for value in column_occupancy
            ),
        }
    )


output_file.parent.mkdir(
    parents=True,
    exist_ok=True
)

with open(output_file, "w", newline="") as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=rows[0].keys(),
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(rows)


print(f"Alignment QC completed: {output_file}")
