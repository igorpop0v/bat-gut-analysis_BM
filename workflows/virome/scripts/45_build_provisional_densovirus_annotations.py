#!/usr/bin/env python3

import csv
import sys
from pathlib import Path
from urllib.parse import quote


def read_fasta(path):
    records = {}
    name = None
    sequence = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(sequence).upper()
                name = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
    if name is not None:
        records[name] = "".join(sequence).upper()
    return records


def read_tsv(path):
    with open(path, newline="") as handle:
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
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


if len(sys.argv) != 5:
    raise SystemExit(
        "Usage: 45_build_provisional_densovirus_annotations.py "
        "CURATED_FASTA ANNOTATION_MANIFEST ORF_REVIEW OUTPUT_DIR"
    )

fasta_path = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
review_path = Path(sys.argv[3])
output_dir = Path(sys.argv[4])
output_dir.mkdir(parents=True, exist_ok=True)

sequences = read_fasta(fasta_path)
manifest = read_tsv(manifest_path)
review = read_tsv(review_path)
metadata = {row["analysis_id"]: row for row in manifest}

gff_path = output_dir / "current_study_densovirus_provisional_annotations.gff3"
tbl_path = output_dir / "current_study_densovirus_provisional_features.tbl"
protein_path = output_dir / "current_study_densovirus_provisional_proteins.faa"
cds_path = output_dir / "current_study_densovirus_provisional_cds.fna"
summary_path = output_dir / "current_study_densovirus_provisional_annotation_summary.tsv"

intact_rows = [row for row in review if row["annotation_status"] == "intact_candidate"]
unresolved_by_genome = {
    row["analysis_id"]: row["unresolved_orfs"]
    for row in manifest
}

feature_counts = {genome_id: 0 for genome_id in sequences}

with (
    open(gff_path, "w") as gff,
    open(tbl_path, "w") as tbl,
    open(protein_path, "w") as proteins,
    open(cds_path, "w") as cds_handle,
):
    gff.write("##gff-version 3\n")
    for genome_id, sequence in sequences.items():
        gff.write(f"##sequence-region {genome_id} 1 {len(sequence)}\n")

    rows_by_genome = {}
    for row in intact_rows:
        rows_by_genome.setdefault(row["analysis_id"], []).append(row)

    for genome_id in sequences:
        tbl.write(f">Feature {genome_id}\n")
        genome_rows = sorted(
            rows_by_genome.get(genome_id, []),
            key=lambda row: int(row["feature_start"]),
        )

        for row in genome_rows:
            start = int(row["feature_start"])
            end = int(row["feature_end"])
            strand = row["strand"]
            orf = row["orf"]
            functional_class = row["functional_class"]
            isolate = metadata[genome_id]["isolate"]

            nucleotide = sequences[genome_id][start - 1:end]
            if strand == "-":
                nucleotide = reverse_complement(nucleotide)

            protein = translate(nucleotide)
            if len(nucleotide) % 3 != 0:
                raise ValueError(f"CDS length is not divisible by three: {genome_id} {orf}")
            if not protein.startswith("M"):
                raise ValueError(f"CDS lacks an ATG start codon: {genome_id} {orf}")
            if not protein.endswith("*"):
                raise ValueError(f"CDS lacks a terminal stop codon: {genome_id} {orf}")
            if "*" in protein[:-1]:
                raise ValueError(f"CDS contains an internal stop codon: {genome_id} {orf}")

            protein = protein[:-1]
            feature_id = f"{genome_id}_{orf}"
            attributes = (
                f"ID={quote(feature_id)};Name={quote(orf)};"
                f"product={quote(orf)};note={quote(functional_class)}"
            )
            gff.write(
                f"{genome_id}\tRefSeq-guided\tCDS\t{start}\t{end}\t.\t"
                f"{strand}\t0\t{attributes}\n"
            )

            if strand == "+":
                table_start, table_end = start, end
            else:
                table_start, table_end = end, start

            tbl.write(f"{table_start}\t{table_end}\tCDS\n")
            tbl.write(f"\t\t\tproduct\t{orf}\n")
            tbl.write(f"\t\t\tnote\t{functional_class}; RefSeq-guided provisional annotation\n")

            proteins.write(
                f">{feature_id} [isolate={isolate}] [product={orf}]\n{wrap(protein)}\n"
            )
            cds_handle.write(
                f">{feature_id} [isolate={isolate}] [product={orf}]\n{wrap(nucleotide)}\n"
            )
            feature_counts[genome_id] += 1

summary_rows = []
for genome_id in sequences:
    unresolved = unresolved_by_genome[genome_id]
    summary_rows.append({
        "analysis_id": genome_id,
        "isolate": metadata[genome_id]["isolate"],
        "sequence_length_nt": len(sequences[genome_id]),
        "intact_cds_written": feature_counts[genome_id],
        "unresolved_orfs_not_annotated": unresolved,
        "annotation_status": (
            "provisional_complete_expected_orf_set"
            if unresolved == "None"
            else "provisional_with_unresolved_orf"
        ),
    })

with open(summary_path, "w", newline="") as handle:
    columns = [
        "analysis_id", "isolate", "sequence_length_nt", "intact_cds_written",
        "unresolved_orfs_not_annotated", "annotation_status",
    ]
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(summary_rows)

print("Provisional densovirus annotations prepared.")
print(f"Sequences: {len(sequences)}")
print(f"Intact CDS written: {sum(feature_counts.values())}")
print(f"GFF3: {gff_path}")
print(f"NCBI feature table draft: {tbl_path}")
print(f"Proteins: {protein_path}")
print(f"CDS nucleotide sequences: {cds_path}")
print(f"Summary: {summary_path}")
