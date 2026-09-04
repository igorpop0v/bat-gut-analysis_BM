#!/usr/bin/env python3

"""Prepare current-study densovirus contigs that require read validation."""

import csv
import sys
from collections import defaultdict
from pathlib import Path


if len(sys.argv) != 6:
    print(
        "Usage: 40_prepare_orf_validation_dataset.py "
        "<genomes.fasta> <genome_manifest.tsv> <orf_audit.tsv> "
        "<output.fasta> <output_manifest.tsv>"
    )
    sys.exit(1)


fasta_file = Path(sys.argv[1])
genome_manifest_file = Path(sys.argv[2])
audit_file = Path(sys.argv[3])
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


sequences = read_fasta(fasta_file)

with genome_manifest_file.open() as handle:
    manifest_rows = list(csv.DictReader(handle, delimiter="\t"))

manifest_by_id = {
    row["analysis_id"]: row
    for row in manifest_rows
}

problems = defaultdict(list)

with audit_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        if (
            row["dataset"] == "current_study"
            and row["annotation_status"] != "intact_candidate"
        ):
            problems[row["genome_id"]].append(
                f'{row["orf"]}:{row["annotation_status"]}'
            )


output_fasta.parent.mkdir(parents=True, exist_ok=True)
output_manifest.parent.mkdir(parents=True, exist_ok=True)

output_rows = []

with output_fasta.open("w") as fasta_handle:
    for genome_id in sorted(problems):
        metadata = manifest_by_id[genome_id]
        original_id = metadata["original_id"]
        sequence = sequences[genome_id]

        # The established read-mapping script obtains the sample ID from
        # the part of this original contig ID preceding the vertical bar.
        fasta_handle.write(f">{original_id}\n")

        for position in range(0, len(sequence), 80):
            fasta_handle.write(sequence[position:position + 80] + "\n")

        output_rows.append(
            {
                "genome_id": genome_id,
                "mapping_reference_id": original_id,
                "sample_id": metadata["sample_id"],
                "group": metadata["group"],
                "length_nt": len(sequence),
                "problem_orfs": ";".join(problems[genome_id])
            }
        )


columns = [
    "genome_id",
    "mapping_reference_id",
    "sample_id",
    "group",
    "length_nt",
    "problem_orfs"
]

with output_manifest.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=columns,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(output_rows)


print("ORF read-validation dataset prepared.")
print(f"Current-study genomes requiring validation: {len(output_rows)}")
print(f"FASTA: {output_fasta}")
print(f"Manifest: {output_manifest}")
