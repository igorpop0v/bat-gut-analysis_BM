#!/usr/bin/env python3

"""Prepare a separate Nr11 sequence with the read-supported nt 1997 deletion."""

import csv
import sys
from pathlib import Path


if len(sys.argv) != 4:
    print(
        "Usage: 42_prepare_nr11_read_supported_candidate.py "
        "<densoviruses.fasta> <genome_manifest.tsv> <output_directory>"
    )
    sys.exit(1)


fasta_file = Path(sys.argv[1])
manifest_file = Path(sys.argv[2])
output_dir = Path(sys.argv[3])
output_dir.mkdir(parents=True, exist_ok=True)

target_id = "Current_Nr11"
deleted_position = 1997
expected_base = "G"


def read_fasta(path):
    records = {}
    sequence_id = None

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                sequence_id = line[1:].split()[0]
                records[sequence_id] = ""
            else:
                records[sequence_id] += line.upper()

    return records


def write_fasta(records, path):
    with path.open("w") as handle:
        for sequence_id, sequence in records.items():
            handle.write(f">{sequence_id}\n")

            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


records = read_fasta(fasta_file)

if target_id not in records:
    raise ValueError(f"Sequence not found: {target_id}")

original_sequence = records[target_id]
observed_base = original_sequence[deleted_position - 1]

if observed_base != expected_base:
    raise ValueError(
        f"Expected {expected_base} at {target_id}:{deleted_position}, "
        f"but found {observed_base}"
    )

corrected_sequence = (
    original_sequence[:deleted_position - 1]
    + original_sequence[deleted_position:]
)

corrected_records = dict(records)
corrected_records[target_id] = corrected_sequence

candidate_fasta = output_dir / "densoviruses_nr11_del1997_candidate.fasta"
candidate_manifest = (
    output_dir / "densoviruses_nr11_del1997_candidate_manifest.tsv"
)
mapping_fasta = output_dir / "nr11_del1997_mapping_reference.fasta"
correction_log = output_dir / "nr11_del1997_correction.tsv"

write_fasta(corrected_records, candidate_fasta)

with manifest_file.open() as handle:
    manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
    manifest_columns = list(manifest_rows[0])

target_manifest_row = None

for row in manifest_rows:
    if row["analysis_id"] == target_id:
        target_manifest_row = row
        row["length_nt"] = str(len(corrected_sequence))
        row["notes"] = (
            row.get("notes", "")
            + "; candidate read-supported deletion of G at position 1997 "
            + "(4/5 Q30, MAPQ20 reads)"
        ).strip("; ")

if target_manifest_row is None:
    raise ValueError(f"Manifest row not found: {target_id}")

with candidate_manifest.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=manifest_columns,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(manifest_rows)

mapping_reference_id = (
    target_manifest_row["original_id"] + "_read_supported_del1997"
)

write_fasta({mapping_reference_id: corrected_sequence}, mapping_fasta)

with correction_log.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(
        [
            "genome_id",
            "operation",
            "original_position",
            "deleted_base",
            "supporting_reads",
            "total_high_quality_reads",
            "support_pct",
            "original_length_nt",
            "candidate_length_nt",
        ]
    )
    writer.writerow(
        [
            target_id,
            "delete",
            deleted_position,
            expected_base,
            4,
            5,
            "80.0",
            len(original_sequence),
            len(corrected_sequence),
        ]
    )


print("Nr11 read-supported candidate prepared.")
print(f"Original length: {len(original_sequence)} nt")
print(f"Candidate length: {len(corrected_sequence)} nt")
print(f"Candidate FASTA: {candidate_fasta}")
print(f"Candidate manifest: {candidate_manifest}")
print(f"Mapping reference: {mapping_fasta}")
print(f"Correction log: {correction_log}")
