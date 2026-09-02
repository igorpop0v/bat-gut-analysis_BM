#!/usr/bin/env python3

import csv
import sys
from pathlib import Path


if len(sys.argv) != 5:
    print(
        "Usage: 23_prepare_densovirus_phylogeny_dataset.py "
        "<external_references.fasta> "
        "<study_genomes.fasta> "
        "<study_manifest.tsv> "
        "<output_directory>"
    )
    sys.exit(1)


external_fasta = Path(sys.argv[1])
study_fasta = Path(sys.argv[2])
study_manifest_file = Path(sys.argv[3])
output_directory = Path(sys.argv[4])

output_directory.mkdir(parents=True, exist_ok=True)

combined_fasta = output_directory / "densovirus_phylogeny_genomes.fasta"
combined_manifest = output_directory / "densovirus_phylogeny_manifest.tsv"


def read_fasta(path):
    records = []
    header = None
    sequence = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append(
                        {
                            "id": header.split()[0],
                            "description": header,
                            "sequence": "".join(sequence).upper()
                        }
                    )

                header = line[1:]
                sequence = []

            else:
                sequence.append(line)

    if header is not None:
        records.append(
            {
                "id": header.split()[0],
                "description": header,
                "sequence": "".join(sequence).upper()
            }
        )

    return records


external_records = read_fasta(external_fasta)
all_study_records = read_fasta(study_fasta)

with study_manifest_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    study_manifest = {
        row["analysis_id"]: row
        for row in reader
    }
    original_fields = reader.fieldnames


allowed_study_datasets = {
    "current_study",
    "previous_study"
}

study_records = [
    record
    for record in all_study_records
    if (
        record["id"] in study_manifest
        and study_manifest[record["id"]]["dataset"]
        in allowed_study_datasets
    )
]


external_ids = {record["id"] for record in external_records}
study_ids = {record["id"] for record in study_records}

duplicate_ids = external_ids.intersection(study_ids)

if duplicate_ids:
    raise ValueError(
        "Duplicate sequence IDs: " +
        ", ".join(sorted(duplicate_ids))
    )


missing_metadata = study_ids.difference(study_manifest)

if missing_metadata:
    raise ValueError(
        "Study sequences missing from the manifest: " +
        ", ".join(sorted(missing_metadata))
    )


with combined_fasta.open("w") as output_handle:
    for record in external_records + study_records:
        output_handle.write(f">{record['id']}\n")

        sequence = record["sequence"]

        for start in range(0, len(sequence), 80):
            output_handle.write(sequence[start:start + 80] + "\n")


manifest_fields = list(original_fields) + ["description"]

manifest_rows = []

for record in external_records:
    accession = record["id"]

    if accession.startswith("NC_"):
        dataset = "refseq_reference"
    else:
        dataset = "genbank_reference"

    manifest_rows.append(
        {
            "analysis_id": accession,
            "dataset": dataset,
            "original_id": accession,
            "sample_id": "NA",
            "group": "NA",
            "status": "complete_reference",
            "length_nt": len(record["sequence"]),
            "source_file": external_fasta.name,
            "notes": "External reference genome",
            "original_orientation": "not_standardized",
            "operation": "downloaded",
            "description": record["description"]
        }
    )


for record in study_records:
    row = study_manifest[record["id"]].copy()
    row["description"] = row["notes"]

    expected_length = int(row["length_nt"])
    observed_length = len(record["sequence"])

    if expected_length != observed_length:
        raise ValueError(
            f"Length mismatch for {record['id']}: "
            f"manifest={expected_length}, FASTA={observed_length}"
        )

    manifest_rows.append(row)


with combined_manifest.open("w", newline="") as output_handle:
    writer = csv.DictWriter(
        output_handle,
        fieldnames=manifest_fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(manifest_rows)


print("Densovirus phylogeny dataset prepared.")
print(f"External reference genomes: {len(external_records)}")
print(f"Previous-study genomes: {sum(r['dataset'] == 'previous_study' for r in manifest_rows)}")
print(f"Current-study genomes: {sum(r['dataset'] == 'current_study' for r in manifest_rows)}")
print(f"Total genomes: {len(manifest_rows)}")
print(f"FASTA: {combined_fasta}")
print(f"Manifest: {combined_manifest}")
