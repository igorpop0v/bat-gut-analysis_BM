#!/usr/bin/env python3

import csv
import sys
from pathlib import Path


if len(sys.argv) != 7:
    print(
        "Usage: 18_prepare_refseq_oriented_dataset.py "
        "<study_sequences.fasta> "
        "<orientation_summary.tsv> "
        "<refseq.fasta> "
        "<study_manifest.tsv> "
        "<output.fasta> "
        "<output_manifest.tsv>"
    )
    sys.exit(1)


study_fasta = Path(sys.argv[1])
orientation_file = Path(sys.argv[2])
refseq_fasta = Path(sys.argv[3])
study_manifest_file = Path(sys.argv[4])
output_fasta = Path(sys.argv[5])
output_manifest = Path(sys.argv[6])


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


study_sequences = read_fasta(study_fasta)
refseq_sequences = read_fasta(refseq_fasta)

if len(refseq_sequences) != 1:
    raise ValueError("The RefSeq FASTA must contain exactly one sequence.")


orientations = {}

with orientation_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        orientations[row["query_id"]] = row["dominant_orientation"]


missing_orientations = sorted(
    set(study_sequences) - set(orientations)
)

if missing_orientations:
    raise ValueError(
        "Missing orientation information for: "
        + ", ".join(missing_orientations)
    )


study_metadata = {}

with study_manifest_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        study_metadata[row["analysis_id"]] = row


output_fasta.parent.mkdir(parents=True, exist_ok=True)
output_manifest.parent.mkdir(parents=True, exist_ok=True)


manifest_columns = [
    "analysis_id",
    "dataset",
    "original_id",
    "sample_id",
    "group",
    "status",
    "length_nt",
    "source_file",
    "notes",
    "original_orientation",
    "operation"
]


with output_fasta.open("w") as fasta_handle, \
        output_manifest.open("w", newline="") as manifest_handle:

    writer = csv.DictWriter(
        manifest_handle,
        fieldnames=manifest_columns,
        delimiter="\t"
    )

    writer.writeheader()

    for sequence_id, sequence in study_sequences.items():
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

        metadata = study_metadata[sequence_id]

        writer.writerow(
            {
                "analysis_id": sequence_id,
                "dataset": metadata["dataset"],
                "original_id": metadata["original_id"],
                "sample_id": metadata["sample_id"],
                "group": metadata["group"],
                "status": metadata["status"],
                "length_nt": len(standardized_sequence),
                "source_file": metadata["source_file"],
                "notes": metadata["notes"],
                "original_orientation": orientation,
                "operation": operation
            }
        )

    refseq_original_id, refseq_sequence = next(
        iter(refseq_sequences.items())
    )

    refseq_analysis_id = "RefSeq_NC_031450.1"

    fasta_handle.write(f">{refseq_analysis_id}\n")

    for position in range(0, len(refseq_sequence), 80):
        fasta_handle.write(
            refseq_sequence[position:position + 80] + "\n"
        )

    writer.writerow(
        {
            "analysis_id": refseq_analysis_id,
            "dataset": "external_reference",
            "original_id": refseq_original_id,
            "sample_id": "NA",
            "group": "NA",
            "status": "public_refseq_reference",
            "length_nt": len(refseq_sequence),
            "source_file": refseq_fasta.name,
            "notes": "Parus major densovirus isolate PmDNV-JL",
            "original_orientation": "Forward",
            "operation": "unchanged"
        }
    )


print(f"Saved FASTA: {output_fasta}")
print(f"Saved manifest: {output_manifest}")
print(f"Study sequences: {len(study_sequences)}")
print("External RefSeq sequences: 1")
print(f"Total sequences: {len(study_sequences) + 1}")
