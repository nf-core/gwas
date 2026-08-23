#!/usr/bin/env python3

import csv
import gzip
import hashlib
import json
import shutil
import sys
from pathlib import Path


REQUIRED_COLUMNS = ["SNPID", "CHR", "POS", "EA", "NEA", "STATUS", "EAF", "BETA", "SE", "P", "N"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path):
    with path.open("rb") as handle:
        is_gzip = handle.read(2) == b"\\x1f\\x8b"
    return gzip.open(path, "rt", encoding="utf-8", newline="") if is_gzip else path.open(
        "rt", encoding="utf-8", newline=""
    )


def inspect_table(path: Path):
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError("canonical summary-statistics table is empty") from error
        repeated = sorted({column for column in header if header.count(column) > 1})
        if repeated:
            raise ValueError(f"canonical summary-statistics header repeats columns: {', '.join(repeated)}")
        missing = [column for column in REQUIRED_COLUMNS if column not in header]
        if missing:
            raise ValueError(f"canonical summary-statistics header is missing required columns: {', '.join(missing)}")
        row_count = 0
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"canonical summary-statistics row {line_number} has {len(row)} fields; expected {len(header)}"
                )
            row_count += 1
        if row_count == 0:
            raise ValueError("canonical summary-statistics table contains no variant rows")
    return header, row_count


def write_canonical(candidate: Path, output: Path):
    with candidate.open("rb") as handle:
        is_gzip = handle.read(2) == b"\\x1f\\x8b"
    if is_gzip:
        shutil.copyfile(candidate, output)
        return "byte_preserving_copy"
    with candidate.open("rb") as source, output.open("wb") as raw_output:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_output, mtime=0) as compressed:
            shutil.copyfileobj(source, compressed)
    return "deterministic_gzip"


def main():
    candidate = Path($candidate_literal)
    source = Path($source_literal)
    output = Path($output_literal)
    provenance_path = Path($provenance_literal)
    metadata = json.loads($metadata_literal)

    header, row_count = inspect_table(candidate)
    serialization = write_canonical(candidate, output)
    provenance = {
        "summary_statistics_id": metadata["summary_statistics_id"],
        "canonical_contract_version": "nfcore_gwas_canonical_v1",
        "trait_id": metadata["trait_id"],
        "trait_type": metadata["trait_type"],
        "genome_build": metadata["build"],
        "ancestry": metadata["ancestry"],
        "source_kind": metadata["source_kind"],
        "source_mode": metadata["source_mode"],
        "source_format": metadata["source_format"],
        "source_method": metadata["source_method"],
        "source_release": metadata.get("source_release"),
        "source_name": metadata["source_name"],
        "source_sha256": sha256(source),
        "producer_analysis_id": metadata.get("producer_analysis_id"),
        "producer_association_method": metadata.get("producer_association_method"),
        "canonical_name": output.name,
        "canonical_sha256": sha256(output),
        "canonical_columns": header,
        "variant_rows": row_count,
        "transformation": metadata["transformation"],
        "serialization": serialization,
        "harmonization": metadata.get("harmonization"),
        "access_constraints": metadata.get("access_constraints"),
    }
    for field in ["population_prevalence", "sample_prevalence"]:
        if metadata.get(field) is not None:
            provenance[field] = metadata[field]
    with provenance_path.open("w", encoding="utf-8") as handle:
        json.dump(provenance, handle, sort_keys=True, indent=2)
        handle.write("\\n")

    with Path("versions.yml").open("w", encoding="utf-8") as versions:
        versions.write('"$task.process":\\n')
        versions.write(f"    python: {sys.version.split()[0]}\\n")


if __name__ == "__main__":
    sys.exit(main())
