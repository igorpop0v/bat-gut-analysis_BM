#!/usr/bin/env python3

import csv
import sys
from pathlib import Path


input_dir = Path(sys.argv[1])
output_dir = Path(sys.argv[2])

output_dir.mkdir(
    parents=True,
    exist_ok=True
)


def read_fasta(fasta_file):

    records = {}
    sequence_id = None
    sequence = []

    with open(fasta_file) as handle:

        for line in handle:

            line = line.strip()

            if line.startswith(">"):

                if sequence_id is not None:
                    records[sequence_id] = "".join(sequence)

                sequence_id = line[1:].split()[0]
                sequence = []

            elif line:
                sequence.append(line)

    if sequence_id is not None:
        records[sequence_id] = "".join(sequence)

    return records


def write_fasta(records, output_file):

    with open(output_file, "w") as handle:

        for sequence_id, sequence in records.items():

            handle.write(f">{sequence_id}\n")

            for position in range(0, len(sequence), 80):
                handle.write(sequence[position:position + 80] + "\n")


summary_rows = []
excluded_rows = []

datasets = {
    "primary": {
        "nsp1": "primary",
        "sp1": "primary",
    },
    "study_inclusive": {
        "nsp1": "sensitivity",
        "sp1": "primary",
    },
    "sensitivity": {
        "nsp1": "sensitivity",
        "sp1": "sensitivity",
    },
}

for dataset, marker_sets in datasets.items():

    nsp1_file = (
        input_dir /
        f"nsp1_{marker_sets['nsp1']}_gt70.faa"
    )

    sp1_file = (
        input_dir /
        f"sp1_{marker_sets['sp1']}_gt70.faa"
    )

    nsp1 = read_fasta(nsp1_file)
    sp1 = read_fasta(sp1_file)

    nsp1_length = len(next(iter(nsp1.values())))
    sp1_length = len(next(iter(sp1.values())))

    shared_ids = sorted(
        set(nsp1) & set(sp1)
    )

    combined = {
        sequence_id: nsp1[sequence_id] + sp1[sequence_id]
        for sequence_id in shared_ids
    }

    dataset_dir = output_dir / dataset

    dataset_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    fasta_output = (
        dataset_dir /
        f"densovirus_{dataset}_nsp1_sp1.faa"
    )

    partition_output = (
        dataset_dir /
        f"densovirus_{dataset}_partitions.nex"
    )

    write_fasta(
        combined,
        fasta_output
    )

    with open(partition_output, "w") as handle:

        handle.write("#nexus\n")
        handle.write("begin sets;\n")
        handle.write(
            f"    charset NSP1 = 1-{nsp1_length};\n"
        )
        handle.write(
            f"    charset SP1 = "
            f"{nsp1_length + 1}-"
            f"{nsp1_length + sp1_length};\n"
        )
        handle.write("end;\n")

    summary_rows.append(
        {
            "dataset": dataset,
            "nsp1_sequences": len(nsp1),
            "sp1_sequences": len(sp1),
            "shared_sequences": len(shared_ids),
            "nsp1_length": nsp1_length,
            "sp1_length": sp1_length,
            "concatenated_length": (
                nsp1_length + sp1_length
            ),
        }
    )

    all_ids = sorted(
        set(nsp1) | set(sp1)
    )

    for sequence_id in all_ids:

        if sequence_id not in shared_ids:

            excluded_rows.append(
                {
                    "dataset": dataset,
                    "sequence_id": sequence_id,
                    "present_in_nsp1": sequence_id in nsp1,
                    "present_in_sp1": sequence_id in sp1,
                    "reason": "Missing one marker",
                }
            )


with open(
    output_dir / "partitioned_alignment_summary.tsv",
    "w",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=summary_rows[0].keys(),
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(summary_rows)


with open(
    output_dir / "excluded_from_partitioned_alignments.tsv",
    "w",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=excluded_rows[0].keys(),
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(excluded_rows)


print("Partitioned alignments prepared.")
print(f"Results: {output_dir}")
