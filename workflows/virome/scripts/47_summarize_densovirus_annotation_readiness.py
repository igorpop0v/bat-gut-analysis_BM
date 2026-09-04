#!/usr/bin/env python3

"""Summarize final annotation QC for current-study densovirus sequences."""

import csv
import sys
from pathlib import Path


if len(sys.argv) != 6:
    raise SystemExit(
        "Usage: 47_summarize_densovirus_annotation_readiness.py "
        "CURATED_FASTA ANNOTATION_MANIFEST ANNOTATION_SUMMARY "
        "PROTEIN_VALIDATION OUTPUT_DIR"
    )


fasta_path = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
summary_path = Path(sys.argv[3])
validation_path = Path(sys.argv[4])
output_dir = Path(sys.argv[5])
output_dir.mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def write_fasta(path, records, selected_ids):
    with path.open("w") as handle:
        for sequence_id in selected_ids:
            sequence = records[sequence_id]
            handle.write(f">{sequence_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


sequences = read_fasta(fasta_path)
manifest = {row["analysis_id"]: row for row in read_tsv(manifest_path)}
annotation_summary = {
    row["analysis_id"]: row for row in read_tsv(summary_path)
}
validation = read_tsv(validation_path)

if set(sequences) != set(manifest) or set(sequences) != set(annotation_summary):
    raise ValueError("FASTA, manifest and annotation summary IDs do not match")

failed_proteins = [
    row for row in validation if row["validation_status"] != "PASS"
]
if failed_proteins:
    raise ValueError(
        f"Protein validation contains {len(failed_proteins)} non-PASS records"
    )

validated_by_genome = {sequence_id: 0 for sequence_id in sequences}
for row in validation:
    genome_id = row["query_id"].rsplit("_ORF", 1)[0]
    if genome_id not in validated_by_genome:
        raise ValueError(f"Unknown genome in protein validation: {genome_id}")
    validated_by_genome[genome_id] += 1

qc_rows = []
full_orf_ids = []
review_ids = []

for sequence_id in sequences:
    metadata = manifest[sequence_id]
    annotation = annotation_summary[sequence_id]
    annotated_cds = int(annotation["intact_cds_written"])
    validated_cds = validated_by_genome[sequence_id]
    unresolved = annotation["unresolved_orfs_not_annotated"]

    if annotated_cds != validated_cds:
        raise ValueError(
            f"Annotated and validated CDS counts differ for {sequence_id}"
        )

    if unresolved == "None" and validated_cds == 5:
        readiness = "complete_expected_orf_set"
        recommendation = "retain_all_five_validated_cds"
        full_orf_ids.append(sequence_id)
    else:
        readiness = "conservative_annotation_with_unresolved_orf"
        recommendation = "retain_only_validated_cds_and_report_unresolved_orf"
        review_ids.append(sequence_id)

    qc_rows.append(
        {
            "analysis_id": sequence_id,
            "isolate": metadata["isolate"],
            "group": metadata["group"],
            "sequence_length_nt": len(sequences[sequence_id]),
            "expected_orfs": 5,
            "annotated_cds": annotated_cds,
            "validated_cds": validated_cds,
            "unresolved_orfs": unresolved,
            "annotation_readiness": readiness,
            "recommended_annotation_strategy": recommendation,
        }
    )

qc_path = output_dir / "current_study_densovirus_annotation_qc.tsv"
with qc_path.open("w", newline="") as handle:
    columns = list(qc_rows[0])
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(qc_rows)

write_fasta(
    output_dir / "densoviruses_complete_expected_orf_set.fasta",
    sequences,
    full_orf_ids,
)
write_fasta(
    output_dir / "densoviruses_unresolved_orf_review_set.fasta",
    sequences,
    review_ids,
)

(output_dir / "complete_expected_orf_set_ids.txt").write_text(
    "\n".join(full_orf_ids) + "\n"
)
(output_dir / "unresolved_orf_review_set_ids.txt").write_text(
    "\n".join(review_ids) + "\n"
)

print("Densovirus annotation-readiness audit completed.")
print(f"Sequences assessed: {len(sequences)}")
print(f"Validated CDS: {len(validation)}")
print(f"Complete expected ORF set: {len(full_orf_ids)}")
print(f"Conservative annotation with unresolved ORF: {len(review_ids)}")
print(f"Summary: {qc_path}")
