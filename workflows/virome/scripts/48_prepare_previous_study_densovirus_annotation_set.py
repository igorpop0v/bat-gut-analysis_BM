#!/usr/bin/env python3

"""Prepare the previous-study densoviruses for uniform annotation."""

import csv
import sys
from pathlib import Path


if len(sys.argv) != 6:
    raise SystemExit(
        "Usage: 48_prepare_previous_study_densovirus_annotation_set.py "
        "COMBINED_FASTA COMBINED_MANIFEST ORF_AUDIT TREE_METADATA OUTPUT_DIR"
    )


fasta_path = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
audit_path = Path(sys.argv[3])
metadata_path = Path(sys.argv[4])
output_dir = Path(sys.argv[5])
output_dir.mkdir(parents=True, exist_ok=True)


def read_fasta(path):
    records = {}
    sequence_id = None
    sequence = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sequence_id is not None:
                    records[sequence_id] = "".join(sequence)
                sequence_id = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)

    if sequence_id is not None:
        records[sequence_id] = "".join(sequence)

    return records


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, columns):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


sequences = read_fasta(fasta_path)
manifest = [
    row for row in read_tsv(manifest_path)
    if row["dataset"] == "previous_study"
]
audit = [
    row for row in read_tsv(audit_path)
    if row["dataset"] == "previous_study"
]
tree_metadata = {row["Name"]: row for row in read_tsv(metadata_path)}

expected_orfs = ["ORF5", "ORF3", "ORF4", "ORF1", "ORF2"]
audit_by_genome = {}
for row in audit:
    audit_by_genome.setdefault(row["genome_id"], []).append(row)

fasta_output = output_dir / "previous_study_densoviruses_curated.fasta"
manifest_output = output_dir / "previous_study_annotation_manifest.tsv"
review_output = output_dir / "previous_study_orf_review.tsv"
unresolved_output = output_dir / "previous_study_unresolved_orfs.tsv"

manifest_rows = []
review_rows = []
unresolved_rows = []

with fasta_output.open("w") as fasta_handle:
    for row in manifest:
        genome_id = row["analysis_id"]
        if genome_id not in sequences:
            raise ValueError(f"Sequence not found: {genome_id}")
        if genome_id not in tree_metadata:
            raise ValueError(f"Tree metadata not found: {genome_id}")

        sequence = sequences[genome_id]
        metadata = tree_metadata[genome_id]
        isolate = metadata["Full.Name"].replace(
            "Densovirinae sp. isolate ", "", 1
        )
        genome_audit = audit_by_genome.get(genome_id, [])
        if len(genome_audit) != 5:
            raise ValueError(
                f"Expected five ORF assessments for {genome_id}, "
                f"found {len(genome_audit)}"
            )

        status_by_orf = {
            item["orf"]: item["annotation_status"]
            for item in genome_audit
        }
        unresolved = [
            orf for orf in expected_orfs
            if status_by_orf.get(orf) != "intact_candidate"
        ]

        fasta_handle.write(f">{genome_id}\n")
        for start in range(0, len(sequence), 80):
            fasta_handle.write(sequence[start:start + 80] + "\n")

        manifest_rows.append(
            {
                "analysis_id": genome_id,
                "isolate": isolate,
                "collection_year": metadata["Year"],
                "country": metadata["Country"],
                "host": metadata["Host"],
                "length_nt": len(sequence),
                "sequence_status": "previous_study_assembled_sequence",
                "intact_orfs": 5 - len(unresolved),
                "unresolved_orfs": ",".join(unresolved) if unresolved else "None",
                "annotation_readiness": (
                    "all_expected_orfs_intact"
                    if not unresolved
                    else "manual_feature_review_required"
                ),
                "sequence_notes": row["notes"],
            }
        )

        for item in genome_audit:
            review_row = {
                "analysis_id": genome_id,
                "isolate": isolate,
                "orf": item["orf"],
                "functional_class": item["functional_class"],
                "reference_protein_id": item["reference_protein_id"],
                "reference_length_aa": item["reference_length_aa"],
                "combined_coverage_pct": item["combined_coverage_pct"],
                "feature_start": item["feature_start"],
                "feature_end": item["feature_end"],
                "strand": item["strand"],
                "translated_length_aa": item["translated_length_aa"],
                "annotation_status": item["annotation_status"],
                "review_note": item["review_note"],
            }
            review_rows.append(review_row)
            if item["annotation_status"] != "intact_candidate":
                unresolved_rows.append(review_row)

manifest_columns = [
    "analysis_id", "isolate", "collection_year", "country", "host",
    "length_nt", "sequence_status", "intact_orfs", "unresolved_orfs",
    "annotation_readiness", "sequence_notes",
]
review_columns = [
    "analysis_id", "isolate", "orf", "functional_class",
    "reference_protein_id", "reference_length_aa", "combined_coverage_pct",
    "feature_start", "feature_end", "strand", "translated_length_aa",
    "annotation_status", "review_note",
]

write_tsv(manifest_output, manifest_rows, manifest_columns)
write_tsv(review_output, review_rows, review_columns)
write_tsv(unresolved_output, unresolved_rows, review_columns)

print("Previous-study densovirus annotation set prepared.")
print(f"Sequences: {len(manifest_rows)}")
print(
    "All expected ORFs intact: "
    + str(sum(row["unresolved_orfs"] == "None" for row in manifest_rows))
)
print(
    "Manual feature review required: "
    + str(sum(row["unresolved_orfs"] != "None" for row in manifest_rows))
)
print(f"FASTA: {fasta_output}")
print(f"Manifest: {manifest_output}")
print(f"ORF review: {review_output}")
print(f"Unresolved ORFs: {unresolved_output}")
