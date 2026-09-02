#!/usr/bin/env python3

import csv
import sys
from collections import defaultdict
from pathlib import Path


if len(sys.argv) != 4:
    raise SystemExit(
        "Usage: 28_reconstruct_study_marker_proteins.py "
        "<tblastn.tsv> <genome_manifest.tsv> <output_dir>"
    )

blast_file = Path(sys.argv[1])
manifest_file = Path(sys.argv[2])
output_dir = Path(sys.argv[3])
output_dir.mkdir(parents=True, exist_ok=True)

marker_by_query = {
    "YP_009310053.1": "NSP1",
    "YP_009310055.1": "SP1",
}


def query_accession(query_id):
    return query_id.split("|")[0]


def write_fasta(path, rows):
    with path.open("w") as handle:
        for row in rows:
            handle.write(f">{row['genome_id']}\n")
            sequence = row["sequence"]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


with manifest_file.open(newline="") as handle:
    manifest_rows = list(csv.DictReader(handle, delimiter="\t"))

study_manifest = {
    row["analysis_id"]: row
    for row in manifest_rows
    if row["dataset"] in {"current_study", "previous_study"}
}

hits = defaultdict(list)

with blast_file.open() as handle:
    reader = csv.reader(handle, delimiter="\t")

    for fields in reader:
        if not fields:
            continue

        query_id = fields[0]
        genome_id = fields[1]
        query_id_clean = query_accession(query_id)

        if query_id_clean not in marker_by_query:
            continue

        if genome_id not in study_manifest:
            continue

        hits[(genome_id, query_id_clean)].append({
            "qstart": int(fields[2]),
            "qend": int(fields[3]),
            "qlen": int(fields[4]),
            "sstart": int(fields[5]),
            "send": int(fields[6]),
            "sframe": int(fields[7]),
            "pident": float(fields[8]),
            "evalue": float(fields[9]),
            "bitscore": float(fields[10]),
            "qseq": fields[11].upper(),
            "sseq": fields[12].upper(),
        })


def reconstruct_protein(hsp_rows):
    qlen = hsp_rows[0]["qlen"]
    reconstructed = [None] * qlen

    for hsp in sorted(
        hsp_rows,
        key=lambda row: row["bitscore"],
        reverse=True,
    ):
        query_position = hsp["qstart"] - 1

        for query_residue, subject_residue in zip(
            hsp["qseq"],
            hsp["sseq"],
        ):
            if query_residue == "-":
                continue

            if (
                reconstructed[query_position] is None
                and subject_residue not in {"-", "*"}
            ):
                reconstructed[query_position] = subject_residue

            query_position += 1

    covered_positions = sum(
        residue is not None for residue in reconstructed
    )
    coverage_pct = 100 * covered_positions / qlen
    sequence = "".join(
        residue if residue is not None else "X"
        for residue in reconstructed
    )

    return sequence, covered_positions, coverage_pct


result_rows = []

for genome_id, metadata in sorted(study_manifest.items()):
    for query_id, marker in marker_by_query.items():
        hsp_rows = hits.get((genome_id, query_id), [])

        if not hsp_rows:
            result_rows.append({
                "genome_id": genome_id,
                "dataset": metadata["dataset"],
                "group": metadata["group"],
                "marker": marker,
                "query_id": query_id,
                "protein_length_aa": 0,
                "hsp_count": 0,
                "covered_positions": 0,
                "coverage_pct": 0,
                "unknown_positions": 0,
                "best_identity_pct": "NA",
                "total_bitscore": 0,
                "inclusion_status": "excluded",
                "sequence": "",
            })
            continue

        sequence, covered_positions, coverage_pct = (
            reconstruct_protein(hsp_rows)
        )

        if coverage_pct >= 90:
            status = "primary"
        else:
            status = "sensitivity_only"

        result_rows.append({
            "genome_id": genome_id,
            "dataset": metadata["dataset"],
            "group": metadata["group"],
            "marker": marker,
            "query_id": query_id,
            "protein_length_aa": len(sequence),
            "hsp_count": len(hsp_rows),
            "covered_positions": covered_positions,
            "coverage_pct": round(coverage_pct, 3),
            "unknown_positions": sequence.count("X"),
            "best_identity_pct": max(
                row["pident"] for row in hsp_rows
            ),
            "total_bitscore": round(
                sum(row["bitscore"] for row in hsp_rows),
                3,
            ),
            "inclusion_status": status,
            "sequence": sequence,
        })


table_fields = [
    "genome_id",
    "dataset",
    "group",
    "marker",
    "query_id",
    "protein_length_aa",
    "hsp_count",
    "covered_positions",
    "coverage_pct",
    "unknown_positions",
    "best_identity_pct",
    "total_bitscore",
    "inclusion_status",
]

summary_file = output_dir / "study_marker_extraction_summary.tsv"

with summary_file.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=table_fields,
        delimiter="\t",
        extrasaction="ignore",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(result_rows)


for marker in ["NSP1", "SP1"]:
    marker_lower = marker.lower()

    primary_rows = [
        row for row in result_rows
        if row["marker"] == marker
        and row["inclusion_status"] == "primary"
    ]

    sensitivity_rows = [
        row for row in result_rows
        if row["marker"] == marker
        and row["inclusion_status"]
        in {"primary", "sensitivity_only"}
    ]

    write_fasta(
        output_dir / f"study_{marker_lower}_primary.faa",
        primary_rows,
    )

    write_fasta(
        output_dir / f"study_{marker_lower}_sensitivity.faa",
        sensitivity_rows,
    )

    primary_count = len(primary_rows)
    sensitivity_only_count = sum(
        row["marker"] == marker
        and row["inclusion_status"] == "sensitivity_only"
        for row in result_rows
    )
    excluded_count = sum(
        row["marker"] == marker
        and row["inclusion_status"] == "excluded"
        for row in result_rows
    )

    print(f"{marker} primary: {primary_count}")
    print(
        f"{marker} sensitivity_only: "
        f"{sensitivity_only_count}"
    )
    print(f"{marker} excluded: {excluded_count}")


print(f"Summary: {summary_file}")
