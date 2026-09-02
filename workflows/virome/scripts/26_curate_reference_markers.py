#!/usr/bin/env python3

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


if len(sys.argv) != 5:
    raise SystemExit(
        "Usage: 26_curate_reference_markers.py "
        "<proteins.faa> <protein_manifest.tsv> <proteinortho.tsv> <output_dir>"
    )

protein_fasta = Path(sys.argv[1])
manifest_file = Path(sys.argv[2])
proteinortho_file = Path(sys.argv[3])
output_dir = Path(sys.argv[4])
output_dir.mkdir(parents=True, exist_ok=True)


def read_fasta(path):
    sequences = {}
    current_id = None
    parts = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(parts)
                current_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(parts)

    return sequences


def write_fasta(path, records):
    with path.open("w") as handle:
        for accession, sequence in records:
            handle.write(f">{accession}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def normalized_text(row):
    return " ".join(
        [row["gene"], row["product"]]
    ).lower().replace("_", "-")


def is_explicit_nsp1(row):
    text = normalized_text(row)
    patterns = [
        r"\bns-?1\b",
        r"\bnon-?structural protein 1\b",
        r"\bnon-?structural protein ns-?1\b",
        r"\bputative ns-?1\b",
    ]
    return any(re.search(pattern, text) for pattern in patterns)


def is_capsid_candidate(row):
    text = normalized_text(row)
    if "nonstructural" in text or "non-structural" in text:
        return False
    patterns = [
        r"\bcapsid\b",
        r"\bstructural protein\b",
        r"\bviral polypeptide vp[0-9]*\b",
        r"\bviral protein(?: 1-4)?\b",
        r"\bvp[0-9]*\b",
        r"\bcap\b",
    ]
    return any(re.search(pattern, text) for pattern in patterns)


def find_orthogroup(marker_id, rows):
    for row in rows:
        members = []
        for cell in row[3:]:
            if cell != "*":
                members.extend(cell.split(","))
        if marker_id in members:
            return set(members)
    raise ValueError(f"Proteinortho group not found for {marker_id}")


sequences = read_fasta(protein_fasta)

with manifest_file.open(newline="") as handle:
    manifest_rows = list(csv.DictReader(handle, delimiter="\t"))

for row in manifest_rows:
    row["length_aa"] = int(row["length_aa"])

by_accession = defaultdict(list)
by_protein = {}
for row in manifest_rows:
    by_accession[row["accession"]].append(row)
    by_protein[row["protein_analysis_id"]] = row

with proteinortho_file.open(newline="") as handle:
    proteinortho_rows = list(csv.reader(handle, delimiter="\t"))

proteinortho_data = proteinortho_rows[1:]
nsp1_group = find_orthogroup("NC_031450.1__CDS002", proteinortho_data)
sp1_group = find_orthogroup("NC_031450.1__CDS004", proteinortho_data)


def choose_marker(accession, marker):
    proteins = by_accession[accession]

    if marker == "NSP1":
        orthogroup_candidates = [
            row for row in proteins
            if row["protein_analysis_id"] in nsp1_group
        ]
        if orthogroup_candidates:
            selected = max(
                orthogroup_candidates,
                key=lambda row: row["length_aa"],
            )
            source = "Proteinortho NSP1 orthogroup"
        else:
            annotation_candidates = [
                row for row in proteins if is_explicit_nsp1(row)
            ]
            if not annotation_candidates:
                return None, "No confident NSP1 candidate"
            selected = max(
                annotation_candidates,
                key=lambda row: row["length_aa"],
            )
            source = "Explicit GenBank NSP1 annotation"

        if selected["length_aa"] >= 400:
            status = "primary"
            reason = "NSP1 candidate is at least 400 aa"
        else:
            status = "sensitivity_only"
            reason = "NSP1 candidate is shorter than 400 aa"

    else:
        annotation_candidates = [
            row for row in proteins if is_capsid_candidate(row)
        ]
        orthogroup_candidates = [
            row for row in proteins
            if row["protein_analysis_id"] in sp1_group
        ]
        if annotation_candidates:
            selected = max(
                annotation_candidates,
                key=lambda row: row["length_aa"],
            )
            source = "GenBank capsid annotation; longest isoform"
        elif orthogroup_candidates:
            selected = max(
                orthogroup_candidates,
                key=lambda row: row["length_aa"],
            )
            source = "Proteinortho capsid orthogroup; longest isoform"
        else:
            return None, "No confident capsid candidate"

        annotation = normalized_text(selected)
        explicitly_incomplete = (
            "truncated" in annotation or "fragment" in annotation
        )

        if selected["length_aa"] >= 350 and not explicitly_incomplete:
            status = "primary"
            reason = "Longest capsid candidate is at least 350 aa"
        else:
            status = "sensitivity_only"
            reason = "Capsid candidate is short or explicitly incomplete"

    return {
        **selected,
        "marker": marker,
        "selection_source": source,
        "inclusion_status": status,
        "selection_reason": reason,
    }, None


selection_rows = []
for accession in sorted(by_accession):
    for marker in ["NSP1", "SP1"]:
        selected, exclusion_reason = choose_marker(accession, marker)
        if selected is not None:
            selection_rows.append(selected)
        else:
            selection_rows.append({
                "accession": accession,
                "marker": marker,
                "protein_analysis_id": "NA",
                "protein_id": "NA",
                "locus_tag": "NA",
                "gene": "NA",
                "product": "NA",
                "location": "NA",
                "length_aa": 0,
                "selection_source": "NA",
                "inclusion_status": "excluded",
                "selection_reason": exclusion_reason,
            })


output_fields = [
    "accession",
    "marker",
    "protein_analysis_id",
    "protein_id",
    "locus_tag",
    "gene",
    "product",
    "location",
    "length_aa",
    "selection_source",
    "inclusion_status",
    "selection_reason",
]

selection_table = output_dir / "reference_marker_selection.tsv"
with selection_table.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=output_fields,
        delimiter="\t",
        extrasaction="ignore",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(selection_rows)


summary_rows = []
for marker in ["NSP1", "SP1"]:
    marker_rows = [row for row in selection_rows if row["marker"] == marker]
    for status in ["primary", "sensitivity_only", "excluded"]:
        summary_rows.append({
            "marker": marker,
            "status": status,
            "genomes": sum(
                row["inclusion_status"] == status for row in marker_rows
            ),
        })

with (output_dir / "reference_marker_summary.tsv").open(
    "w", newline=""
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["marker", "status", "genomes"],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(summary_rows)


for marker in ["NSP1", "SP1"]:
    marker_lower = marker.lower()
    primary_rows = [
        row for row in selection_rows
        if row["marker"] == marker
        and row["inclusion_status"] == "primary"
    ]
    sensitivity_rows = [
        row for row in selection_rows
        if row["marker"] == marker
        and row["inclusion_status"] in {"primary", "sensitivity_only"}
    ]

    for label, rows in [
        ("primary", primary_rows),
        ("sensitivity", sensitivity_rows),
    ]:
        records = []
        for row in rows:
            protein_id = row["protein_analysis_id"]
            if protein_id not in sequences:
                raise ValueError(f"Sequence missing from FASTA: {protein_id}")
            records.append((row["accession"], sequences[protein_id]))

        write_fasta(
            output_dir / f"external_{marker_lower}_{label}.faa",
            records,
        )


print("Reference marker curation completed.")
print(f"Selection table: {selection_table}")
for row in summary_rows:
    print(f"{row['marker']} {row['status']}: {row['genomes']}")
