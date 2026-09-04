#!/usr/bin/env python3

"""Check read support around split TBLASTN hits in densovirus ORFs."""

import csv
import os
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


if len(sys.argv) != 6:
    print(
        "Usage: 41_validate_orf_junctions.py "
        "<genome_manifest.tsv> <orf_audit.tsv> <tblastn_hits.tsv> "
        "<read_mapping_directory> <output_directory>"
    )
    sys.exit(1)


manifest_file = Path(sys.argv[1]).resolve()
audit_file = Path(sys.argv[2]).resolve()
hits_file = Path(sys.argv[3]).resolve()
mapping_dir = Path(sys.argv[4]).resolve()
output_dir = Path(sys.argv[5]).resolve()
pileup_dir = output_dir / "pileup"

output_dir.mkdir(parents=True, exist_ok=True)
pileup_dir.mkdir(parents=True, exist_ok=True)

docker_image = "bat-gut-densovirus-mapping:1.0"
window_padding = 50


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_id(sequence_id):
    return sequence_id.replace("|", "__")


def parse_pileup_bases(bases):
    """Count mismatches and indel events in a samtools mpileup base string."""

    reference_matches = 0
    mismatches = 0
    insertions = 0
    deletions = 0
    index = 0

    while index < len(bases):
        symbol = bases[index]

        if symbol == "^":
            index += 2
            continue

        if symbol == "$":
            index += 1
            continue

        if symbol in ".,":
            reference_matches += 1
            index += 1
            continue

        if symbol in "ACGTNacgtn":
            mismatches += 1
            index += 1
            continue

        if symbol in "+-":
            is_insertion = symbol == "+"
            index += 1
            length_start = index

            while index < len(bases) and bases[index].isdigit():
                index += 1

            if length_start == index:
                continue

            indel_length = int(bases[length_start:index])
            index += indel_length

            if is_insertion:
                insertions += 1
            else:
                deletions += 1

            continue

        # Deletion placeholders, reference skips and any uncommon symbols.
        index += 1

    return reference_matches, mismatches, insertions, deletions


manifest_rows = read_tsv(manifest_file)
manifest = {row["analysis_id"]: row for row in manifest_rows}

problem_rows = []

for row in read_tsv(audit_file):
    if (
        row["dataset"] == "current_study"
        and row["annotation_status"] != "intact_candidate"
    ):
        problem_rows.append(row)


# TBLASTN columns:
# qseqid sseqid pident length qstart qend qlen
# sstart send slen evalue bitscore sframe
hits = defaultdict(list)

with hits_file.open() as handle:
    reader = csv.reader(handle, delimiter="\t")

    for fields in reader:
        if len(fields) < 13:
            continue

        protein_id = fields[0].split("|")[0]
        genome_id = fields[1]
        hits[(genome_id, protein_id)].append(
            {
                "query_start": int(fields[4]),
                "query_end": int(fields[5]),
                "subject_start": int(fields[7]),
                "subject_end": int(fields[8]),
                "bitscore": float(fields[11]),
                "frame": int(fields[12]),
            }
        )


summary_rows = []
alternative_rows = []

for problem in problem_rows:
    genome_id = problem["genome_id"]
    protein_id = problem["reference_protein_id"]
    metadata = manifest[genome_id]
    mapping_reference_id = metadata["original_id"]
    sample_id = metadata["sample_id"]
    sequence_length = int(metadata["length_nt"])
    file_id = safe_id(mapping_reference_id)

    reference_file = mapping_dir / "references" / f"{file_id}.fasta"
    bam_file = mapping_dir / "bam" / f"{file_id}.bam"

    if not reference_file.exists():
        raise FileNotFoundError(f"Missing mapping reference: {reference_file}")

    if not bam_file.exists():
        raise FileNotFoundError(f"Missing BAM file: {bam_file}")

    protein_hits = hits.get((genome_id, protein_id), [])
    strand = problem["dominant_strand"]
    dominant_hits = [
        hit
        for hit in protein_hits
        if ("+" if hit["frame"] > 0 else "-") == strand
    ]

    if len(dominant_hits) < 2:
        raise ValueError(
            f"Expected at least two dominant-strand HSPs for "
            f"{genome_id} {problem['orf']}"
        )

    dominant_hits.sort(key=lambda hit: hit["query_start"])
    boundary_positions = []

    for left_hit, right_hit in zip(dominant_hits, dominant_hits[1:]):
        boundary_positions.extend(
            [left_hit["subject_end"], right_hit["subject_start"]]
        )

    window_start = max(1, min(boundary_positions) - window_padding)
    window_end = min(sequence_length, max(boundary_positions) + window_padding)
    region = f"{mapping_reference_id}:{window_start}-{window_end}"

    reference_in_container = f"/mapping/references/{reference_file.name}"
    bam_in_container = f"/mapping/bam/{bam_file.name}"

    subprocess.run(
        [
            "docker", "run", "--rm",
            "--user", f"{os.getuid()}:{os.getgid()}",
            "--volume", f"{mapping_dir}:/mapping",
            docker_image,
            "samtools", "faidx", reference_in_container,
        ],
        check=True,
    )

    pileup = subprocess.run(
        [
            "docker", "run", "--rm",
            "--user", f"{os.getuid()}:{os.getgid()}",
            "--volume", f"{mapping_dir}:/mapping",
            docker_image,
            "samtools", "mpileup",
            "-A", "-B", "-q", "20", "-Q", "30", "-d", "10000",
            "-f", reference_in_container,
            "-r", region,
            bam_in_container,
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout

    pileup_file = pileup_dir / f"{genome_id}_{problem['orf']}.pileup.tsv"
    pileup_file.write_text(pileup)

    depths = []
    mismatch_percentages = []
    indel_percentages = []

    for line in pileup.splitlines():
        fields = line.split("\t")

        if len(fields) < 5:
            continue

        position = int(fields[1])
        reference_base = fields[2]
        depth = int(fields[3])
        bases = fields[4]
        matches, mismatches, insertions, deletions = parse_pileup_bases(bases)
        indels = insertions + deletions
        mismatch_pct = 100 * mismatches / depth if depth else 0
        indel_pct = 100 * indels / depth if depth else 0

        depths.append(depth)
        mismatch_percentages.append(mismatch_pct)
        indel_percentages.append(indel_pct)

        if mismatch_pct >= 1 or indel_pct >= 1:
            alternative_rows.append(
                {
                    "genome_id": genome_id,
                    "sample_id": sample_id,
                    "orf": problem["orf"],
                    "position": position,
                    "reference_base": reference_base,
                    "depth": depth,
                    "reference_matches": matches,
                    "mismatches": mismatches,
                    "insertions": insertions,
                    "deletions": deletions,
                    "mismatch_pct": f"{mismatch_pct:.3f}",
                    "indel_pct": f"{indel_pct:.3f}",
                }
            )

    if not depths:
        raise ValueError(f"No pileup positions returned for {genome_id} {problem['orf']}")

    maximum_alternative_pct = max(
        max(mismatch_percentages),
        max(indel_percentages),
    )

    if min(depths) < 20:
        validation_status = "insufficient_local_depth"
    elif maximum_alternative_pct >= 20:
        validation_status = "alternative_allele_requires_review"
    elif maximum_alternative_pct >= 5:
        validation_status = "minor_alternative_signal"
    else:
        validation_status = "assembled_junction_supported"

    summary_rows.append(
        {
            "genome_id": genome_id,
            "sample_id": sample_id,
            "group": metadata["group"],
            "orf": problem["orf"],
            "annotation_status": problem["annotation_status"],
            "frames": problem["frames"],
            "junction_coordinates": ",".join(map(str, boundary_positions)),
            "validation_window": f"{window_start}-{window_end}",
            "positions_examined": len(depths),
            "minimum_depth": min(depths),
            "median_depth": f"{statistics.median(depths):.1f}",
            "maximum_mismatch_pct": f"{max(mismatch_percentages):.3f}",
            "maximum_indel_pct": f"{max(indel_percentages):.3f}",
            "validation_status": validation_status,
        }
    )


summary_columns = [
    "genome_id",
    "sample_id",
    "group",
    "orf",
    "annotation_status",
    "frames",
    "junction_coordinates",
    "validation_window",
    "positions_examined",
    "minimum_depth",
    "median_depth",
    "maximum_mismatch_pct",
    "maximum_indel_pct",
    "validation_status",
]

alternative_columns = [
    "genome_id",
    "sample_id",
    "orf",
    "position",
    "reference_base",
    "depth",
    "reference_matches",
    "mismatches",
    "insertions",
    "deletions",
    "mismatch_pct",
    "indel_pct",
]

summary_file = output_dir / "orf_junction_read_validation.tsv"
alternative_file = output_dir / "orf_junction_alternative_alleles.tsv"

with summary_file.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=summary_columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(summary_rows)

with alternative_file.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=alternative_columns,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(alternative_rows)


print("ORF-junction read validation completed.")
print(f"ORFs examined: {len(summary_rows)}")
print(f"Summary: {summary_file}")
print(f"Alternative alleles >= 1%: {alternative_file}")
