#!/usr/bin/env python3

"""Audit RefSeq-guided densovirus ORFs before GenBank annotation."""

import csv
import sys
from collections import defaultdict
from pathlib import Path


if len(sys.argv) != 6:
    print(
        "Usage: 39_audit_densovirus_orfs.py "
        "<genomes.fasta> <genome_manifest.tsv> <refseq_proteins.tsv> "
        "<tblastn_hits.tsv> <output_directory>"
    )
    sys.exit(1)


fasta_file = Path(sys.argv[1])
genome_manifest_file = Path(sys.argv[2])
protein_manifest_file = Path(sys.argv[3])
hits_file = Path(sys.argv[4])
output_dir = Path(sys.argv[5])
output_dir.mkdir(parents=True, exist_ok=True)


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
                    records[sequence_id] = "".join(sequence_parts).upper()

                sequence_id = line[1:].split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line)

    if sequence_id is not None:
        records[sequence_id] = "".join(sequence_parts).upper()

    return records


def reverse_complement(sequence):
    table = str.maketrans("ACGTN", "TGCAN")
    return sequence.translate(table)[::-1]


codons = [
    a + b + c
    for a in "TCAG"
    for b in "TCAG"
    for c in "TCAG"
]

amino_acids = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)

codon_table = dict(zip(codons, amino_acids))


def translate(sequence):
    return "".join(
        codon_table.get(sequence[i:i + 3], "X")
        for i in range(0, len(sequence) - 2, 3)
    )


sequences = read_fasta(fasta_file)

with genome_manifest_file.open() as handle:
    genome_rows = list(csv.DictReader(handle, delimiter="\t"))

genome_metadata = {
    row["analysis_id"]: row
    for row in genome_rows
}

with protein_manifest_file.open() as handle:
    protein_rows = list(csv.DictReader(handle, delimiter="\t"))

protein_metadata = {
    row["protein_id"]: row
    for row in protein_rows
}


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
                "identity": float(fields[2]),
                "alignment_length": int(fields[3]),
                "query_start": int(fields[4]),
                "query_end": int(fields[5]),
                "query_length": int(fields[6]),
                "subject_start": int(fields[7]),
                "subject_end": int(fields[8]),
                "subject_length": int(fields[9]),
                "e_value": float(fields[10]),
                "bitscore": float(fields[11]),
                "frame": int(fields[12])
            }
        )


output_rows = []

for genome_id in sequences:
    if genome_id not in genome_metadata:
        raise ValueError(f"No manifest row for {genome_id}")

    genome = sequences[genome_id]
    metadata = genome_metadata[genome_id]

    for protein_id, protein in protein_metadata.items():
        protein_hits = hits.get((genome_id, protein_id), [])

        row = {
            "genome_id": genome_id,
            "dataset": metadata.get("dataset", "NA"),
            "group": metadata.get("group", "NA"),
            "reference_protein_id": protein_id,
            "analysis_alias": protein["analysis_alias"],
            "orf": protein["product"],
            "functional_class": protein["note"],
            "reference_length_aa": protein["length_aa"],
            "hsp_count": len(protein_hits),
            "best_hsp_coverage_pct": "NA",
            "combined_coverage_pct": "NA",
            "dominant_strand": "NA",
            "frames": "NA",
            "feature_start": "NA",
            "feature_end": "NA",
            "strand": "NA",
            "start_codon": "NA",
            "stop_codon": "NA",
            "translated_length_aa": "NA",
            "internal_stop_count": "NA",
            "annotation_status": "absent",
            "review_note": "No significant TBLASTN hit"
        }

        if not protein_hits:
            output_rows.append(row)
            continue

        best_hit = max(protein_hits, key=lambda item: item["bitscore"])
        query_length = best_hit["query_length"]
        best_coverage = (
            100
            * (best_hit["query_end"] - best_hit["query_start"] + 1)
            / query_length
        )

        strand_scores = defaultdict(float)

        for hit in protein_hits:
            strand = "+" if hit["frame"] > 0 else "-"
            strand_scores[strand] += hit["bitscore"]

        dominant_strand = max(strand_scores, key=strand_scores.get)
        dominant_hits = [
            hit
            for hit in protein_hits
            if ("+" if hit["frame"] > 0 else "-") == dominant_strand
        ]

        covered_positions = set()

        for hit in dominant_hits:
            covered_positions.update(
                range(hit["query_start"], hit["query_end"] + 1)
            )

        combined_coverage = 100 * len(covered_positions) / query_length
        frames = sorted({hit["frame"] for hit in dominant_hits})

        row.update(
            {
                "best_hsp_coverage_pct": f"{best_coverage:.3f}",
                "combined_coverage_pct": f"{combined_coverage:.3f}",
                "dominant_strand": dominant_strand,
                "frames": ",".join(map(str, frames))
            }
        )

        # A single HSP spanning the complete reference protein permits an
        # exact provisional CDS call. The three terminal nucleotides are
        # added because TBLASTN does not align the stop codon.
        if best_coverage >= 95 and len(frames) == 1:
            start = best_hit["subject_start"]
            end = best_hit["subject_end"]
            expected_coding_nt = query_length * 3
            observed_coding_nt = abs(start - end) + 1

            if observed_coding_nt != expected_coding_nt:
                row["annotation_status"] = "manual_review"
                row["review_note"] = (
                    "Complete protein hit, but genomic span is inconsistent "
                    "with one uninterrupted CDS"
                )
                output_rows.append(row)
                continue

            if dominant_strand == "+":
                feature_start = start
                feature_end = end + 3
                cds = genome[feature_start - 1:feature_end]
            else:
                feature_start = end - 3
                feature_end = start
                cds = reverse_complement(
                    genome[feature_start - 1:feature_end]
                )

            if feature_start < 1 or feature_end > len(genome):
                row["annotation_status"] = "manual_review"
                row["review_note"] = (
                    "Complete protein hit reaches a sequence boundary; "
                    "the terminal stop codon cannot be verified"
                )
                output_rows.append(row)
                continue

            protein_translation = translate(cds)
            internal_stops = protein_translation[:-1].count("*")
            start_codon = cds[:3]
            stop_codon = cds[-3:]
            terminal_stop = protein_translation.endswith("*")

            row.update(
                {
                    "feature_start": feature_start,
                    "feature_end": feature_end,
                    "strand": dominant_strand,
                    "start_codon": start_codon,
                    "stop_codon": stop_codon,
                    "translated_length_aa": len(protein_translation) - 1,
                    "internal_stop_count": internal_stops
                }
            )

            if (
                start_codon == "ATG"
                and terminal_stop
                and internal_stops == 0
            ):
                row["annotation_status"] = "intact_candidate"
                row["review_note"] = (
                    "Complete single-frame protein match with verified "
                    "ATG start and terminal stop codon"
                )
            else:
                row["annotation_status"] = "manual_review"
                row["review_note"] = (
                    "Provisional CDS does not pass start/stop/internal-stop "
                    "checks"
                )

        elif combined_coverage >= 95 and len(frames) > 1:
            row["annotation_status"] = "possible_frameshift"
            row["review_note"] = (
                "Protein is covered by HSPs in more than one reading frame"
            )
        elif combined_coverage >= 95:
            row["annotation_status"] = "split_alignment_manual_review"
            row["review_note"] = (
                "Protein is covered by multiple HSPs; inspect the junction"
            )
        else:
            row["annotation_status"] = "partial_manual_review"
            row["review_note"] = (
                "Protein match is incomplete relative to the RefSeq protein"
            )

        output_rows.append(row)


output_columns = [
    "genome_id",
    "dataset",
    "group",
    "reference_protein_id",
    "analysis_alias",
    "orf",
    "functional_class",
    "reference_length_aa",
    "hsp_count",
    "best_hsp_coverage_pct",
    "combined_coverage_pct",
    "dominant_strand",
    "frames",
    "feature_start",
    "feature_end",
    "strand",
    "start_codon",
    "stop_codon",
    "translated_length_aa",
    "internal_stop_count",
    "annotation_status",
    "review_note"
]

audit_file = output_dir / "densovirus_orf_integrity_audit.tsv"

with audit_file.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=output_columns,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(output_rows)


status_counts = defaultdict(int)

for row in output_rows:
    status_counts[row["annotation_status"]] += 1

summary_file = output_dir / "densovirus_orf_integrity_summary.tsv"

with summary_file.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["annotation_status", "features"])

    for status in sorted(status_counts):
        writer.writerow([status, status_counts[status]])


print("Densovirus ORF integrity audit completed.")
print(f"Genomes: {len(sequences)}")
print(f"Protein profiles per genome: {len(protein_metadata)}")
print(f"Feature assessments: {len(output_rows)}")
print(f"Audit: {audit_file}")
print(f"Summary: {summary_file}")
