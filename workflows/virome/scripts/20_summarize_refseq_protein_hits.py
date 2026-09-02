#!/usr/bin/env python3

import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path


if len(sys.argv) != 6:
    print(
        "Usage: 20_summarize_refseq_protein_hits.py "
        "<tblastn.tsv> <protein_manifest.tsv> "
        "<genome_manifest.tsv> <long_output.tsv> "
        "<summary_output.tsv>"
    )
    sys.exit(1)


blast_file = Path(sys.argv[1])
protein_manifest_file = Path(sys.argv[2])
genome_manifest_file = Path(sys.argv[3])
long_output = Path(sys.argv[4])
summary_output = Path(sys.argv[5])


def interval_coverage(intervals):
    covered_positions = set()

    for start, end in intervals:
        first = min(start, end)
        last = max(start, end)

        covered_positions.update(range(first, last + 1))

    return len(covered_positions)


proteins = []

with protein_manifest_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        proteins.append(row)


genomes = []

with genome_manifest_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        genomes.append(row)


hits = defaultdict(list)

with blast_file.open() as handle:
    reader = csv.reader(handle, delimiter="\t")

    for row in reader:
        query_id = row[0].split("|")[0]

        subject_id = row[1]
        subject_parts = subject_id.split("|")

        if (
            len(subject_parts) >= 3
            and subject_parts[0] in {"ref", "gb", "emb", "dbj"}
        ):
            genome_id = subject_parts[1]
        else:
            genome_id = subject_id

        hits[(query_id, genome_id)].append(
            {
                "identity_pct": float(row[2]),
                "alignment_length": int(row[3]),
                "query_start": int(row[4]),
                "query_end": int(row[5]),
                "query_length": int(row[6]),
                "subject_start": int(row[7]),
                "subject_end": int(row[8]),
                "subject_length": int(row[9]),
                "e_value": float(row[10]),
                "bitscore": float(row[11]),
                "subject_frame": int(row[12])
            }
        )

long_rows = []

for protein in proteins:
    protein_id = protein["protein_id"]
    protein_length = int(protein["length_aa"])

    for genome in genomes:
        genome_id = genome["analysis_id"]
        current_hits = hits.get((protein_id, genome_id), [])

        if current_hits:
            covered_aa = interval_coverage(
                [
                    (
                        hit["query_start"],
                        hit["query_end"]
                    )
                    for hit in current_hits
                ]
            )

            covered_nt = interval_coverage(
                [
                    (
                        hit["subject_start"],
                        hit["subject_end"]
                    )
                    for hit in current_hits
                ]
            )

            alignment_sum = sum(
                hit["alignment_length"]
                for hit in current_hits
            )

            weighted_identity = sum(
                hit["identity_pct"] *
                hit["alignment_length"]
                for hit in current_hits
            ) / alignment_sum

            positive_bitscore = sum(
                hit["bitscore"]
                for hit in current_hits
                if hit["subject_frame"] > 0
            )

            negative_bitscore = sum(
                hit["bitscore"]
                for hit in current_hits
                if hit["subject_frame"] < 0
            )

            protein_coverage = (
                100 * covered_aa / protein_length
            )

            if protein_coverage >= 95:
                hit_status = "complete_like"
            elif protein_coverage >= 80:
                hit_status = "near_complete"
            else:
                hit_status = "partial"

            dominant_strand = (
                "Positive"
                if positive_bitscore >= negative_bitscore
                else "Negative"
            )

            row = {
                "protein_id": protein_id,
                "analysis_alias": protein["analysis_alias"],
                "genome_id": genome_id,
                "dataset": genome["dataset"],
                "group": genome["group"],
                "protein_length_aa": protein_length,
                "hsp_count": len(current_hits),
                "best_identity_pct": round(
                    max(
                        hit["identity_pct"]
                        for hit in current_hits
                    ),
                    3
                ),
                "weighted_identity_pct": round(
                    weighted_identity,
                    3
                ),
                "protein_covered_aa": covered_aa,
                "protein_coverage_pct": round(
                    protein_coverage,
                    3
                ),
                "genome_covered_nt": covered_nt,
                "positive_frame_bitscore": round(
                    positive_bitscore,
                    1
                ),
                "negative_frame_bitscore": round(
                    negative_bitscore,
                    1
                ),
                "dominant_subject_strand": dominant_strand,
                "minimum_e_value": min(
                    hit["e_value"]
                    for hit in current_hits
                ),
                "total_bitscore": round(
                    sum(
                        hit["bitscore"]
                        for hit in current_hits
                    ),
                    1
                ),
                "hit_status": hit_status
            }

        else:
            row = {
                "protein_id": protein_id,
                "analysis_alias": protein["analysis_alias"],
                "genome_id": genome_id,
                "dataset": genome["dataset"],
                "group": genome["group"],
                "protein_length_aa": protein_length,
                "hsp_count": 0,
                "best_identity_pct": "NA",
                "weighted_identity_pct": "NA",
                "protein_covered_aa": 0,
                "protein_coverage_pct": 0,
                "genome_covered_nt": 0,
                "positive_frame_bitscore": 0,
                "negative_frame_bitscore": 0,
                "dominant_subject_strand": "NA",
                "minimum_e_value": "NA",
                "total_bitscore": 0,
                "hit_status": "absent"
            }

        long_rows.append(row)


long_columns = [
    "protein_id",
    "analysis_alias",
    "genome_id",
    "dataset",
    "group",
    "protein_length_aa",
    "hsp_count",
    "best_identity_pct",
    "weighted_identity_pct",
    "protein_covered_aa",
    "protein_coverage_pct",
    "genome_covered_nt",
    "positive_frame_bitscore",
    "negative_frame_bitscore",
    "dominant_subject_strand",
    "minimum_e_value",
    "total_bitscore",
    "hit_status"
]


long_output.parent.mkdir(parents=True, exist_ok=True)

with long_output.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=long_columns,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(long_rows)


summary_rows = []

for protein in proteins:
    protein_rows = [
        row
        for row in long_rows
        if row["protein_id"] == protein["protein_id"]
    ]

    coverage_values = [
        float(row["protein_coverage_pct"])
        for row in protein_rows
    ]

    summary_rows.append(
        {
            "protein_id": protein["protein_id"],
            "analysis_alias": protein["analysis_alias"],
            "length_aa": protein["length_aa"],
            "genomes_total": len(protein_rows),
            "complete_like_hits": sum(
                row["hit_status"] == "complete_like"
                for row in protein_rows
            ),
            "near_complete_hits": sum(
                row["hit_status"] == "near_complete"
                for row in protein_rows
            ),
            "partial_hits": sum(
                row["hit_status"] == "partial"
                for row in protein_rows
            ),
            "absent_hits": sum(
                row["hit_status"] == "absent"
                for row in protein_rows
            ),
            "minimum_coverage_pct": round(
                min(coverage_values),
                3
            ),
            "median_coverage_pct": round(
                statistics.median(coverage_values),
                3
            ),
            "maximum_coverage_pct": round(
                max(coverage_values),
                3
            )
        }
    )


summary_columns = [
    "protein_id",
    "analysis_alias",
    "length_aa",
    "genomes_total",
    "complete_like_hits",
    "near_complete_hits",
    "partial_hits",
    "absent_hits",
    "minimum_coverage_pct",
    "median_coverage_pct",
    "maximum_coverage_pct"
]


with summary_output.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=summary_columns,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(summary_rows)


print(f"Detailed results: {long_output}")
print(f"Marker summary: {summary_output}")
