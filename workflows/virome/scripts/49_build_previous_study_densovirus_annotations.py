#!/usr/bin/env python3

"""Build conservative RefSeq-guided annotations for previous-study sequences."""

import csv
import sys
from pathlib import Path
from urllib.parse import quote


if len(sys.argv) != 5:
    raise SystemExit(
        "Usage: 49_build_previous_study_densovirus_annotations.py "
        "CURATED_FASTA ANNOTATION_MANIFEST ORF_REVIEW OUTPUT_DIR"
    )


fasta_path = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
review_path = Path(sys.argv[3])
output_dir = Path(sys.argv[4])
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
                    records[sequence_id] = "".join(sequence).upper()
                sequence_id = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
    if sequence_id is not None:
        records[sequence_id] = "".join(sequence).upper()
    return records


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


codons = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
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


def wrap(sequence, width=80):
    return "\n".join(
        sequence[start:start + width]
        for start in range(0, len(sequence), width)
    )


sequences = read_fasta(fasta_path)
manifest = read_tsv(manifest_path)
review = read_tsv(review_path)
metadata = {row["analysis_id"]: row for row in manifest}

if set(sequences) != set(metadata):
    raise ValueError("FASTA and annotation manifest IDs do not match")

intact_rows = [
    row for row in review
    if row["annotation_status"] == "intact_candidate"
]
rows_by_genome = {}
for row in intact_rows:
    rows_by_genome.setdefault(row["analysis_id"], []).append(row)

gff_path = output_dir / "previous_study_densovirus_provisional_annotations.gff3"
tbl_path = output_dir / "previous_study_densovirus_provisional_features.tbl"
protein_path = output_dir / "previous_study_densovirus_provisional_proteins.faa"
cds_path = output_dir / "previous_study_densovirus_provisional_cds.fna"
summary_path = output_dir / "previous_study_densovirus_provisional_annotation_summary.tsv"

feature_counts = {sequence_id: 0 for sequence_id in sequences}

with (
    gff_path.open("w") as gff,
    tbl_path.open("w") as tbl,
    protein_path.open("w") as proteins,
    cds_path.open("w") as cds_handle,
):
    gff.write("##gff-version 3\n")
    for sequence_id, sequence in sequences.items():
        gff.write(f"##sequence-region {sequence_id} 1 {len(sequence)}\n")

    for sequence_id, sequence in sequences.items():
        tbl.write(f">Feature {sequence_id}\n")
        genome_rows = sorted(
            rows_by_genome.get(sequence_id, []),
            key=lambda row: int(row["feature_start"]),
        )

        for row in genome_rows:
            start = int(row["feature_start"])
            end = int(row["feature_end"])
            strand = row["strand"]
            orf = row["orf"]
            functional_class = row["functional_class"]
            isolate = metadata[sequence_id]["isolate"]

            nucleotide = sequence[start - 1:end]
            if strand == "-":
                nucleotide = reverse_complement(nucleotide)

            protein = translate(nucleotide)
            if len(nucleotide) % 3 != 0:
                raise ValueError(
                    f"CDS length is not divisible by three: {sequence_id} {orf}"
                )
            if not protein.startswith("M"):
                raise ValueError(f"CDS lacks ATG start: {sequence_id} {orf}")
            if not protein.endswith("*"):
                raise ValueError(f"CDS lacks terminal stop: {sequence_id} {orf}")
            if "*" in protein[:-1]:
                raise ValueError(f"CDS has internal stop: {sequence_id} {orf}")

            protein = protein[:-1]
            feature_id = f"{sequence_id}_{orf}"
            attributes = (
                f"ID={quote(feature_id)};Name={quote(orf)};"
                f"product={quote(orf)};note={quote(functional_class)}"
            )
            gff.write(
                f"{sequence_id}\tRefSeq-guided\tCDS\t{start}\t{end}\t.\t"
                f"{strand}\t0\t{attributes}\n"
            )

            table_start, table_end = (
                (start, end) if strand == "+" else (end, start)
            )
            tbl.write(f"{table_start}\t{table_end}\tCDS\n")
            tbl.write(f"\t\t\tproduct\t{orf}\n")
            tbl.write(
                f"\t\t\tnote\t{functional_class}; "
                "RefSeq-guided provisional annotation\n"
            )

            proteins.write(
                f">{feature_id} [isolate={isolate}] [product={orf}]\n"
                f"{wrap(protein)}\n"
            )
            cds_handle.write(
                f">{feature_id} [isolate={isolate}] [product={orf}]\n"
                f"{wrap(nucleotide)}\n"
            )
            feature_counts[sequence_id] += 1

summary_rows = []
for sequence_id in sequences:
    unresolved = metadata[sequence_id]["unresolved_orfs"]
    summary_rows.append(
        {
            "analysis_id": sequence_id,
            "isolate": metadata[sequence_id]["isolate"],
            "sequence_length_nt": len(sequences[sequence_id]),
            "intact_cds_written": feature_counts[sequence_id],
            "unresolved_orfs_not_annotated": unresolved,
            "annotation_status": (
                "provisional_complete_expected_orf_set"
                if unresolved == "None"
                else "provisional_with_unresolved_orf"
            ),
        }
    )

with summary_path.open("w", newline="") as handle:
    columns = list(summary_rows[0])
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(summary_rows)

print("Previous-study provisional annotations prepared.")
print(f"Sequences: {len(sequences)}")
print(f"Intact CDS written: {sum(feature_counts.values())}")
print(f"GFF3: {gff_path}")
print(f"NCBI feature table draft: {tbl_path}")
print(f"Proteins: {protein_path}")
print(f"CDS nucleotide sequences: {cds_path}")
print(f"Summary: {summary_path}")
