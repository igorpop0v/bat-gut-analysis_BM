#!/usr/bin/env python3

"""Validate the read-supported Nr11 deletion after remapping."""

import csv
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path


if len(sys.argv) != 3:
    print(
        "Usage: 43_validate_nr11_corrected_consensus.py "
        "<corrected_read_mapping_directory> <output_directory>"
    )
    sys.exit(1)


mapping_dir = Path(sys.argv[1]).resolve()
output_dir = Path(sys.argv[2]).resolve()
output_dir.mkdir(parents=True, exist_ok=True)

image = "bat-gut-densovirus-mapping:1.0"
sequence_id = "RNA_S12046Nr11|k119_2287_read_supported_del1997"
file_id = sequence_id.replace("|", "__")
reference_file = mapping_dir / "references" / f"{file_id}.fasta"
bam_file = mapping_dir / "bam" / f"{file_id}.bam"
region = f"{sequence_id}:1990-2002"
anchor_position = 1996

if not reference_file.exists():
    raise FileNotFoundError(f"Missing reference FASTA: {reference_file}")

if not bam_file.exists():
    raise FileNotFoundError(f"Missing BAM file: {bam_file}")

reference_in_container = f"/mapping/references/{reference_file.name}"
bam_in_container = f"/mapping/bam/{bam_file.name}"

subprocess.run(
    [
        "docker", "run", "--rm",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "--volume", f"{mapping_dir}:/mapping",
        image,
        "samtools", "faidx", reference_in_container,
    ],
    check=True,
)

pileup = subprocess.run(
    [
        "docker", "run", "--rm",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "--volume", f"{mapping_dir}:/mapping",
        image,
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

raw_file = output_dir / "nr11_corrected_region.pileup.tsv"
raw_file.write_text(pileup)


def parse_bases(bases):
    matches = 0
    mismatches = 0
    insertions = Counter()
    deletions = Counter()
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
            matches += 1
            index += 1
            continue

        if symbol in "ACGTNacgtn":
            mismatches += 1
            index += 1
            continue

        if symbol in "+-":
            operation = symbol
            index += 1
            length_start = index

            while index < len(bases) and bases[index].isdigit():
                index += 1

            length = int(bases[length_start:index])
            allele = bases[index:index + length].upper()
            index += length

            if operation == "+":
                insertions[allele] += 1
            else:
                deletions[allele] += 1

            continue

        index += 1

    return matches, mismatches, insertions, deletions


position_rows = []
anchor_row = None

for line in pileup.splitlines():
    fields = line.split("\t")

    if len(fields) < 5:
        continue

    position = int(fields[1])
    depth = int(fields[3])
    matches, mismatches, insertions, deletions = parse_bases(fields[4])

    row = {
        "position": position,
        "reference_base": fields[2],
        "depth": depth,
        "reference_matches": matches,
        "mismatches": mismatches,
        "insertions": sum(insertions.values()),
        "insertion_alleles": ";".join(
            f"{allele}:{count}" for allele, count in sorted(insertions.items())
        ) or "NA",
        "deletions": sum(deletions.values()),
        "deletion_alleles": ";".join(
            f"{allele}:{count}" for allele, count in sorted(deletions.items())
        ) or "NA",
    }
    position_rows.append(row)

    if position == anchor_position:
        anchor_row = row

if anchor_row is None:
    raise ValueError(f"No pileup result at anchor position {anchor_position}")

original_g_reads = 0

for item in anchor_row["insertion_alleles"].split(";"):
    if item.startswith("G:"):
        original_g_reads += int(item.split(":")[1])

corrected_reads = anchor_row["depth"] - anchor_row["insertions"]
corrected_support_pct = (
    100 * corrected_reads / anchor_row["depth"]
    if anchor_row["depth"]
    else 0
)

if corrected_reads >= 3 and corrected_support_pct >= 75:
    decision = "corrected_consensus_supported"
else:
    decision = "ambiguous_consensus"

position_file = output_dir / "nr11_corrected_region_summary.tsv"
decision_file = output_dir / "nr11_corrected_consensus_validation.tsv"

with position_file.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(position_rows[0]),
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(position_rows)

with decision_file.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(
        [
            "sequence_id",
            "anchor_position",
            "high_quality_depth",
            "corrected_no_insertion_reads",
            "original_G_insertion_reads",
            "corrected_support_pct",
            "validation_status",
        ]
    )
    writer.writerow(
        [
            sequence_id,
            anchor_position,
            anchor_row["depth"],
            corrected_reads,
            original_g_reads,
            f"{corrected_support_pct:.1f}",
            decision,
        ]
    )

print("Nr11 corrected consensus validation completed.")
print(f"Decision: {decision}")
print(f"Summary: {decision_file}")
print(f"Regional evidence: {position_file}")
