#!/usr/bin/env python3

import csv
import gzip
import hashlib
import json
import math
import sys
from pathlib import Path


INPUT = Path(${input_literal})
OUTPUT = Path(${output_literal})
PREPARATION = Path(${preparation_literal})
META = json.loads(${metadata_literal})
PROCESS = ${task_process_literal}
REQUIRED = ["SNPID", "EA", "NEA", "EAF", "BETA", "SE", "N"]


def open_text(path):
    with path.open("rb") as handle:
        magic = handle.read(2)
    return gzip.open(path, "rt", encoding="utf-8", newline="") if magic == bytes([31, 139]) else path.open("r", encoding="utf-8", newline="")


def finite_number(value, field, row_number):
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"row {row_number} field {field!r} is not numeric: {value!r}") from error
    if not math.isfinite(number):
        raise ValueError(f"row {row_number} field {field!r} is not finite: {value!r}")
    return number


row_count = 0
ambiguous_count = 0
with open_text(INPUT) as source, OUTPUT.open("w", encoding="utf-8", newline="") as destination:
    reader = csv.DictReader(source, delimiter="\\t")
    missing = [column for column in REQUIRED if column not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"canonical summary statistics are missing required columns: {', '.join(missing)}")
    writer = csv.writer(destination, delimiter="\\t", lineterminator="\\n")
    writer.writerow(["Predictor", "A1", "A2", "Z", "n", "A1Freq"])
    for row_number, row in enumerate(reader, start=2):
        beta = finite_number(row["BETA"], "BETA", row_number)
        se = finite_number(row["SE"], "SE", row_number)
        n = finite_number(row["N"], "N", row_number)
        eaf = finite_number(row["EAF"], "EAF", row_number)
        if se <= 0:
            raise ValueError(f"row {row_number} field 'SE' must be greater than zero")
        if n <= 0:
            raise ValueError(f"row {row_number} field 'N' must be greater than zero")
        if not 0 <= eaf <= 1:
            raise ValueError(f"row {row_number} field 'EAF' must be within [0,1]")
        allele_pair = {row["EA"].upper(), row["NEA"].upper()}
        if allele_pair in [{"A", "T"}, {"C", "G"}]:
            ambiguous_count += 1
        writer.writerow([row["SNPID"], row["EA"], row["NEA"], format(beta / se, ".17g"), row["N"], row["EAF"]])
        row_count += 1

if row_count == 0:
    raise ValueError("canonical summary statistics contain no data rows")

sha256 = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
preparation = {
    "schema_version": "1.0",
    "summary_statistics_id": META["summary_statistics_id"],
    "canonical_contract_version": META.get("canonical_contract_version"),
    "conversion": "canonical_beta_se_to_ldak_z",
    "rows_written": row_count,
    "strand_ambiguous_rows_present": ambiguous_count,
    "strand_ambiguous_policy_owner": "LDAK native --allow-ambiguous default",
    "output_sha256": sha256,
    "execution_owner": PROCESS,
}
PREPARATION.write_text(json.dumps(preparation, indent=2, sort_keys=True) + "\\n", encoding="utf-8")
Path("versions.yml").write_text(f'"{PROCESS}":\\n    python: {sys.version.split()[0]}\\n', encoding="utf-8")
