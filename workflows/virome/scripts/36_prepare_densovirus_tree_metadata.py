#!/usr/bin/env python3

"""Prepare Host, Country, collection Year and Study annotations for trees."""

import csv
import re
import sys
from pathlib import Path


def qualifier(record, name):
    match = re.search(
        rf'/{re.escape(name)}="(.*?)"',
        record,
        flags=re.DOTALL,
    )
    if match is None:
        return "ND"
    return re.sub(r"\s+", " ", match.group(1)).strip()


def definition(record):
    match = re.search(
        r"^DEFINITION\s+(.*?)(?=^ACCESSION\s+)",
        record,
        flags=re.MULTILINE | re.DOTALL,
    )
    if match is None:
        return "ND"
    value = re.sub(r"\s+", " ", match.group(1)).strip()
    return value.rstrip(".")


def infer_host(full_name):
    name = re.sub(r"^\S+\s+", "", full_name)

    if re.search(r"bat[- ]associated", name, flags=re.IGNORECASE):
        return "bat"
    if re.search(r"mosquito", name, flags=re.IGNORECASE):
        return "Mosquito"

    host = re.split(
        r"\s+(?:denso\S*|ambi\S*|itera\S*|metallo\S*|hepan\S*)",
        name,
        maxsplit=1,
        flags=re.IGNORECASE,
    )[0]

    if host == name or host.startswith("MAG:"):
        return "ND"

    return host.strip()


def icon_query(host):
    if host in {"ND", ""}:
        return "ND"

    rules = [
        (r"bat|Nyctalus|Miniopterus", "Eptesicus"),
        (r"Anas", "Anas"),
        (r"Homo sapiens", "Homo sapiens"),
        (r"Culex|Aedes|Anopheles|Mosquito", "Culex"),
        (r"Penaeus|Fenneropenaeus", "Penaeus monodon"),
        (r"Cherax", "Decapoda"),
        (r"Actinonaias", "Unionidae"),
        # The current PhyloPic API has no direct entry for Bactericera
        # trigonica, but it does contain the parent superfamily Psylloidea.
        (r"Bactericera", "Psylloidea"),
        (r"Diaphorina", "Diaphorina citri"),
        (r"Asterias", "Asterias"),
        (r"Acheta|cricket", "Acheta domestica"),
        (r"Myzus", "Aphididae"),
        (r"Sibine|Casphalia", "Acharia stimulea"),
        (r"Dendrolimus|Erinnyis", "Sphinx ligustri"),
        (r"Mythimna|Helicoverpa", "Catocala nupta"),
        (r"Galleria", "Eriocephala calthella"),
        (r"Junonia", "Vanessa atalanta"),
        (r"Periplaneta", "Periplaneta americana"),
        (r"Bombus", "Bombus"),
        (r"Dendroctonus", "Dendroctonus"),
        (r"Parus", "Parus"),
        (r"Perigonia", "Dilophonotini"),
        (r"Pseudoplusia", "Spodoptera frugiperda"),
        (r"Neodiprion", "Xyela julii"),
        (r"Diatraea", "Crambus praefectellus"),
        (r"Solenopsis", "Solenopsis invicta"),
        (r"Zophobas|Tenebrio", "Tribolium castaneum"),
        (r"Bombyx", "Bombyx"),
        (r"Danaus", "Danaus plexippus"),
        (r"Blattella", "Blattella germanica"),
        (r"Planococcus", "Planococcus citri"),
    ]

    for pattern, replacement in rules:
        if re.search(pattern, host, flags=re.IGNORECASE):
            return replacement

    return host


def country_code(country):
    country_root = country.split(":", 1)[0].strip()
    mapping = {
        "USA": "us",
        "Russia": "ru",
        "China": "cn",
        "France": "fr",
        "Brazil": "br",
        "India": "in",
        "Canada": "ca",
        "Australia": "au",
        "Belgium": "be",
        "Tanzania": "tz",
        "Israel": "il",
        "Argentina": "ar",
        "Colombia": "co",
        "Croatia": "hr",
        "Taiwan": "tw",
    }
    return mapping.get(country_root, "ND")


def read_genbank_metadata(path):
    text = Path(path).read_text(encoding="utf-8")
    records = {}

    for record in re.split(r"\n//\s*\n", text):
        version_match = re.search(r"^VERSION\s+(\S+)", record, re.MULTILINE)
        if version_match is None:
            continue

        accession = version_match.group(1)
        full_name = f"{accession} {definition(record)}"
        host = qualifier(record, "host")
        country = qualifier(record, "geo_loc_name")
        if country == "ND":
            country = qualifier(record, "country")

        date = qualifier(record, "collection_date")
        year_match = re.search(r"(?:19|20)\d{2}", date)
        year = year_match.group(0) if year_match else "ND"

        if host == "ND":
            host = infer_host(full_name)

        records[accession] = {
            "Full.Name": full_name,
            "Country": country,
            "Year": year,
            "Host": host,
            "Host.icon.query": icon_query(host),
            "flag_code": country_code(country),
        }

    return records


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    if len(sys.argv) not in {5, 6}:
        raise SystemExit(
            "Usage: 36_prepare_densovirus_tree_metadata.py "
            "REFERENCE_RECORDS.gb PHYLOGENY_MANIFEST.tsv "
            "CUSTOM_LABELS.tsv OUTPUT.tsv [CURRENT_STUDY_YEAR]"
        )

    genbank_file, manifest_file, labels_file, output_file = sys.argv[1:5]
    # Kept only for compatibility with the earlier command-line interface.
    # Collection years are now assigned explicitly from the study design.
    legacy_current_year = sys.argv[5] if len(sys.argv) == 6 else None

    current_study_years = {
        "ABH": "2023",
        "H": "2024",
        "AAH": "2024",
    }

    focal_reference_studies = {
        "NC_031450.1": "Yang et al., 2016",
        "MW628494.1": "Armién et al., 2023",
    }

    genbank = read_genbank_metadata(genbank_file)
    manifest = read_tsv(manifest_file)
    labels = {
        row["analysis_id"]: row["display_label"]
        for row in read_tsv(labels_file)
    }

    output_rows = []

    for row in manifest:
        analysis_id = row["analysis_id"]
        dataset = row["dataset"]

        if dataset in {"refseq_reference", "genbank_reference"}:
            values = genbank.get(analysis_id)
            if values is None:
                raise ValueError(f"No GenBank metadata for {analysis_id}")

            display_label = re.sub(
                r", complete.*$",
                "",
                values["Full.Name"],
            )
            display_label = re.sub(
                r" densovirus non-structural protein.*$",
                " densovirus",
                display_label,
            )
            group = "ND"
            study = focal_reference_studies.get(analysis_id, "ND")

        else:
            display_label = labels.get(analysis_id, row["description"])
            group = row["group"] if dataset == "current_study" else "ND"

            if dataset == "previous_study":
                collection_year = "2022"
                study = "Popov et al., 2025"
            elif dataset == "current_study":
                if group not in current_study_years:
                    raise ValueError(
                        "No collection year specified for current-study "
                        f"group {group!r}"
                    )
                collection_year = current_study_years[group]
                study = f"This study: {group}"
            else:
                raise ValueError(f"Unexpected study dataset: {dataset}")

            values = {
                "Country": "Russia",
                "Host": "Nyctalus noctula",
                "Host.icon.query": "Eptesicus",
                "flag_code": "ru",
                "Year": collection_year,
            }

        output_rows.append(
            {
                "Name": analysis_id,
                "Full.Name": display_label,
                "Country": values["Country"],
                "Year": values["Year"],
                "Host": values["Host"],
                "Host.icon.query": values["Host.icon.query"],
                "flag_code": values["flag_code"],
                "Dataset": dataset,
                "Group": group,
                "Study": study,
            }
        )

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "Name",
        "Full.Name",
        "Country",
        "Year",
        "Host",
        "Host.icon.query",
        "flag_code",
        "Dataset",
        "Group",
        "Study",
    ]

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    print("Densovirus tree metadata prepared.")
    print(f"Sequences: {len(output_rows)}")
    print(f"Output: {output_path}")
    if legacy_current_year is not None:
        print(
            "Note: the legacy CURRENT_STUDY_YEAR argument was ignored; "
            "collection years were assigned from the experimental design."
        )


if __name__ == "__main__":
    main()
