#!/usr/bin/env python3

import csv
import sys
from pathlib import Path


ASSOCIATIONS = [Path(path) for path in ${associations_literal}]
PROCESS = ${task_process_literal}

# The GWASLab column mapping registered for `ldak_kvik` binds neither of these: `MAF` is a minor-allele
# frequency that the joined `EAF` replaces, and `Wald_Stat` is the statistic whose `Wald_P` is bound instead.
UNBOUND_COLUMNS = ("Wald_Stat", "MAF")
KEY_COLUMNS = ("Predictor", "A1", "A2")


def read_summary(path):
    with path.open(encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\\t")
        return {
            tuple(row[column] for column in KEY_COLUMNS): (row["A1Freq"], row["n"])
            for row in reader
        }


for association in ASSOCIATIONS:
    frequency_and_size = read_summary(association.with_suffix(".summaries"))
    with association.open(encoding="utf-8", newline="") as source, Path(
        association.name.removesuffix(".assoc") + ".harmonisation.tsv"
    ).open("w", encoding="utf-8", newline="") as destination:
        reader = csv.DictReader(source, delimiter="\\t")
        carried = [column for column in reader.fieldnames if column not in UNBOUND_COLUMNS]
        writer = csv.writer(destination, delimiter="\\t", lineterminator="\\n")
        writer.writerow(carried + ["EAF", "N"])
        for row in reader:
            key = tuple(row[column] for column in KEY_COLUMNS)
            effect_allele_frequency, sample_size = frequency_and_size.get(key, ("", ""))
            writer.writerow([row[column] for column in carried] + [effect_allele_frequency, sample_size])

Path("versions.yml").write_text(f'"{PROCESS}":\\n    python: {sys.version.split()[0]}\\n', encoding="utf-8")
