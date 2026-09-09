#!/usr/bin/env python3
"""Convert aligned GWASLab study columns to native MR-MEGA inputs."""

import csv
import gzip
import json
import math
import os
import sys


Z_95 = 1.959963984540054


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", newline="")


def chromosome_code(value):
    text = str(value).removeprefix("chr").removeprefix("CHR")
    return {"X": "23", "Y": "24", "XY": "25", "M": "26", "MT": "26"}.get(text.upper(), text)


def count_text(value):
    number = float(value)
    return str(int(number)) if number.is_integer() else repr(number)


def main():
    aligned_studies = $aligned_studies_literal
    study_names = $study_names_literal
    trait_type = $trait_type_literal
    prefix = $prefix_literal
    binary = trait_type == "binary"

    directory = "{}.mrmega_inputs".format(prefix)
    os.makedirs(directory, exist_ok=True)
    width = max(2, len(str(len(study_names))))
    study_paths = [
        os.path.join(directory, "{}_{}.txt.gz".format(str(index + 1).zfill(width), name))
        for index, name in enumerate(study_names)
    ]
    handles = [gzip.GzipFile(filename="", mode="wb", fileobj=open(path, "wb"), mtime=0) for path in study_paths]
    header = (
        ["MARKERNAME", "EA", "NEA", "EAF", "OR", "OR_95L", "OR_95U", "N", "CHROMOSOME", "POSITION"]
        if binary
        else ["MARKERNAME", "EA", "NEA", "EAF", "BETA", "SE", "N", "CHROMOSOME", "POSITION"]
    )

    try:
        for handle in handles:
            handle.write(("\\t".join(header) + "\\n").encode("utf-8"))
        with open_text(aligned_studies) as source:
            for row in csv.DictReader(source, delimiter="\t"):
                for index, handle in enumerate(handles, start=1):
                    beta = row["BETA_{}".format(index)]
                    error = row["SE_{}".format(index)]
                    frequency = row["EAF_{}".format(index)]
                    size = row["N_{}".format(index)]
                    if "NA" in {beta, error, frequency, size}:
                        continue
                    effect = float(beta)
                    standard_error = float(error)
                    effect_fields = (
                        [
                            repr(math.exp(effect)),
                            repr(math.exp(effect - Z_95 * standard_error)),
                            repr(math.exp(effect + Z_95 * standard_error)),
                        ]
                        if binary
                        else [repr(effect), repr(standard_error)]
                    )
                    fields = [row["META_VARIANT_KEY"], row["EA"], row["NEA"], repr(float(frequency))]
                    fields.extend(effect_fields)
                    fields.extend([count_text(size), chromosome_code(row["CHR"]), row["POS"]])
                    handle.write(("\\t".join(fields) + "\\n").encode("utf-8"))
    finally:
        for handle in handles:
            handle.close()

    with open("{}.filelist.txt".format(prefix), "w", encoding="utf-8", newline="\\n") as handle:
        for path in study_paths:
            handle.write(path + "\\n")
    with open("versions.yml", "w", encoding="utf-8", newline="\\n") as handle:
        handle.write('"$task.process":\\n')
        handle.write("    python: {}\\n".format(sys.version.split()[0]))


if __name__ == "__main__":
    sys.exit(main())
