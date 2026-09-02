#!/usr/bin/env python3

import csv
import sys
from pathlib import Path

from Bio import SeqIO


if len(sys.argv) != 4:
    print(
        "Usage: 24_extract_reference_cds.py "
        "<input.gb> <output.faa> <output.tsv>"
    )
    sys.exit(1)


genbank_file = Path(sys.argv[1])
protein_fasta = Path(sys.argv[2])
protein_manifest = Path(sys.argv[3])


def clean_text(value):
    return " ".join(str(value).replace("\t", " ").split())


protein_fasta.parent.mkdir(parents=True, exist_ok=True)
protein_manifest.parent.mkdir(parents=True, exist_ok=True)

manifest_rows = []
record_count = 0
records_with_proteins = set()
protein_count = 0


with protein_fasta.open("w") as fasta_handle:
    for record in SeqIO.parse(genbank_file, "genbank"):
        record_count += 1
        accession = record.id
        cds_index = 0

        for feature in record.features:
            if feature.type != "CDS":
                continue

            translation = feature.qualifiers.get(
                "translation",
                [None]
            )[0]

            if translation is None:
                continue

            cds_index += 1
            protein_count += 1
            records_with_proteins.add(accession)

            analysis_id = (
                f"{accession}__CDS{cds_index:03d}"
            )

            protein_id = feature.qualifiers.get(
                "protein_id",
                ["NA"]
            )[0]

            locus_tag = feature.qualifiers.get(
                "locus_tag",
                ["NA"]
            )[0]

            gene = feature.qualifiers.get(
                "gene",
                ["NA"]
            )[0]

            product = feature.qualifiers.get(
                "product",
                ["unannotated protein"]
            )[0]

            translation = clean_text(translation)

            fasta_handle.write(f">{analysis_id}\n")

            for start in range(0, len(translation), 80):
                fasta_handle.write(
                    translation[start:start + 80] + "\n"
                )

            manifest_rows.append(
                {
                    "protein_analysis_id": analysis_id,
                    "accession": accession,
                    "protein_id": clean_text(protein_id),
                    "locus_tag": clean_text(locus_tag),
                    "gene": clean_text(gene),
                    "product": clean_text(product),
                    "location": str(feature.location),
                    "length_aa": len(translation)
                }
            )


columns = [
    "protein_analysis_id",
    "accession",
    "protein_id",
    "locus_tag",
    "gene",
    "product",
    "location",
    "length_aa"
]


with protein_manifest.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=columns,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(manifest_rows)


print("Reference CDS extraction completed.")
print(f"GenBank records: {record_count}")
print(
    "Records containing translated CDS: "
    f"{len(records_with_proteins)}"
)
print(f"Translated proteins: {protein_count}")
print(f"Protein FASTA: {protein_fasta}")
print(f"Protein manifest: {protein_manifest}")
