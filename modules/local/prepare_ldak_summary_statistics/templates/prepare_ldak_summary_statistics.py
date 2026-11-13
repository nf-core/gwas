#!/usr/bin/env python3

import csv
import gzip
import sys
from pathlib import Path


INPUT = Path(${input_literal})
OUTPUT = Path(${output_literal})
PROCESS = ${task_process_literal}


with gzip.open(INPUT, "rt", encoding="utf-8", newline="") as source, OUTPUT.open(
    "w", encoding="utf-8", newline=""
) as destination:
    reader = csv.DictReader(source, delimiter="\t")
    writer = csv.writer(destination, delimiter="\t", lineterminator="\\n")
    writer.writerow(["Predictor", "A1", "A2", "Z", "n", "A1Freq"])
    for row in reader:
        writer.writerow(
            [
                row["SNPID"],
                row["EA"],
                row["NEA"],
                float(row["BETA"]) / float(row["SE"]),
                row["N"],
                row["EAF"],
            ]
        )

Path("versions.yml").write_text(f'"{PROCESS}":\\n    python: {sys.version.split()[0]}\\n', encoding="utf-8")
