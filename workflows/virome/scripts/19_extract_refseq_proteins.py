#!/usr/bin/env python3

import csv
import re
import sys
from pathlib import Path


if len(sys.argv) != 4:
    print(
        "Usage: 19_extract_refseq_proteins.py "
        "<input.gb> <output.faa> <output_manifest.tsv>"
    )
    sys.exit(1)


genbank_file = Path(sys.argv[1])
output_fasta = Path(sys.argv[2])
output_manifest = Path(sys.argv[3])


analysis_aliases = {
    "YP_009310052.1": "ORF5_nonstructural",
    "YP_009310053.1": "NSP1_candidate_ORF3",
    "YP_009310054.1": "ORF4_nonstructural",
    "YP_009310055.1": "SP1_candidate_ORF1",
    "YP_009310056.1": "SP2_candidate_ORF2"
}


def extract_qualifier(block, qualifier):
    pattern = rf'/{qualifier}="(.*?)"'
    match = re.search(pattern, block, flags=re.DOTALL)

    if match is None:
        return ""

    return re.sub(
        r"\s+",
        " ",
        match.group(1)
    ).strip()


genbank_text = genbank_file.read_text()

features_text = genbank_text.split(
    "FEATURES", 1
)[1].split(
    "ORIGIN", 1
)[0]

feature_blocks = re.split(
    r"(?=^     \S)",
    features_text,
    flags=re.MULTILINE
)

proteins = []

for block in feature_blocks:
    first_line = block.splitlines()[0]

    if not first_line.startswith("     CDS"):
        continue

    location = first_line[21:].strip()
    protein_id = extract_qualifier(block, "protein_id")
    locus_tag = extract_qualifier(block, "locus_tag")
    product = extract_qualifier(block, "product")
    note = extract_qualifier(block, "note")
    translation = extract_qualifier(block, "translation")
    translation = re.sub(r"\s+", "", translation)

    if not protein_id or not translation:
        continue

    proteins.append(
        {
            "protein_id": protein_id,
            "analysis_alias": analysis_aliases.get(
                protein_id,
                "unassigned"
            ),
            "locus_tag": locus_tag,
            "product": product,
            "note": note,
            "location": location,
            "length_aa": len(translation),
            "translation": translation
        }
    )


output_fasta.parent.mkdir(parents=True, exist_ok=True)
output_manifest.parent.mkdir(parents=True, exist_ok=True)


with output_fasta.open("w") as handle:
    for protein in proteins:
        handle.write(
            f">{protein['protein_id']}|"
            f"{protein['analysis_alias']}\n"
        )

        sequence = protein["translation"]

        for position in range(0, len(sequence), 80):
            handle.write(sequence[position:position + 80] + "\n")


manifest_columns = [
    "protein_id",
    "analysis_alias",
    "locus_tag",
    "product",
    "note",
    "location",
    "length_aa"
]


with output_manifest.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=manifest_columns,
        delimiter="\t"
    )

    writer.writeheader()

    for protein in proteins:
        writer.writerow(
            {
                column: protein[column]
                for column in manifest_columns
            }
        )


print(f"Proteins extracted: {len(proteins)}")
print(f"Protein FASTA: {output_fasta}")
print(f"Protein manifest: {output_manifest}")
