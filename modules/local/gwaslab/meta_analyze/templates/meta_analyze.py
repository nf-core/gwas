#!/usr/bin/env python3
"""Align GWASLab parents and emit native meta-analysis and adapter products."""

import gzip
import json
import os
import sys

# Importing GWASLab 4.1.9 as a non-root container user raises RuntimeError when
# Numba cache=True targets the read-only package directory. Retire when the bioconda::gwaslab pin moves past the fix.
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.abspath(".numba_cache"))
os.environ.setdefault("MPLCONFIGDIR", os.path.abspath(".matplotlib"))

import gwaslab as gl
import pandas as pd
from gwaslab.g_SumstatsMulti import SumstatsMulti


BUILDS = {"GRCh37": "19", "GRCh38": "38"}


def exact(value):
    if pd.isna(value):
        return "NA"
    return repr(float(value))


def write_result(frame, path):
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        na_rep="NA",
        compression={"method": "gzip", "mtime": 0},
    )


def main():
    parents = json.loads(r'''$parents_literal''')
    study_names = json.loads(r'''$study_names_literal''')
    input_format = json.loads(r'''$input_format_literal''')
    genome_build = json.loads(r'''$genome_build_literal''')
    random_effects = json.loads(r'''$random_effects_literal''')
    prefix = json.loads(r'''$prefix_literal''')

    objects = [
        gl.Sumstats(
            path,
            fmt=input_format,
            build=BUILDS[genome_build],
            species="homo sapiens",
            study=study_name,
            verbose=False,
        )
        for path, study_name in zip(parents, study_names)
    ]
    multi = SumstatsMulti(
        objects,
        group_name=prefix,
        build=BUILDS[genome_build],
        engine="pandas",
        merge_by_id=False,
        keep_all_variants=True,
        verbose=False,
    )
    aligned = multi.data.copy()

    fixed = multi.run_meta_analysis(random_effects=False)
    write_result(fixed.data, "{}.fixed.tsv.gz".format(prefix))
    if random_effects:
        random = multi.run_meta_analysis(random_effects=True)
        write_result(random.data, "{}.random.tsv.gz".format(prefix))

    study_count = len(study_names)
    key = (
        aligned["CHR"].astype("string")
        + ":"
        + aligned["POS"].astype("Int64").astype("string")
        + ":"
        + aligned["EA"].astype("string")
        + ":"
        + aligned["NEA"].astype("string")
    )
    aligned.insert(0, "META_VARIANT_KEY", key)

    with open("{}.metasoft.input.txt".format(prefix), "w", encoding="utf-8", newline="\\n") as handle:
        for _, row in aligned.iterrows():
            fields = [str(row["META_VARIANT_KEY"])]
            for index in range(1, study_count + 1):
                beta = row["BETA_{}".format(index)]
                error = row["SE_{}".format(index)]
                if pd.isna(beta) or pd.isna(error):
                    fields.extend(["NA", "NA"])
                else:
                    fields.extend([exact(beta), exact(error)])
            handle.write(" ".join(fields) + "\\n")

    view_columns = ["META_VARIANT_KEY", "SNPID", "CHR", "POS", "EA", "NEA"]
    for index in range(1, study_count + 1):
        view_columns.extend(
            [
                "EAF_{}".format(index),
                "BETA_{}".format(index),
                "SE_{}".format(index),
                "N_{}".format(index),
            ]
        )
    with gzip.GzipFile(
        filename="", mode="wb", fileobj=open("{}.mrmega.tsv.gz".format(prefix), "wb"), mtime=0
    ) as handle:
        handle.write(("\\t".join(view_columns) + "\\n").encode("utf-8"))
        for _, row in aligned.iterrows():
            fields = []
            for column in view_columns:
                value = row[column]
                if column in {"META_VARIANT_KEY", "SNPID", "CHR", "EA", "NEA"}:
                    fields.append("NA" if pd.isna(value) else str(value))
                elif column == "POS":
                    fields.append("NA" if pd.isna(value) else str(int(value)))
                else:
                    fields.append(exact(value))
            handle.write(("\\t".join(fields) + "\\n").encode("utf-8"))

    with open("{}.gwaslab.log".format(prefix), "w", encoding="utf-8", newline="\\n") as handle:
        handle.write(multi.log.log_text)
        handle.write("\\n")

    with open("versions.yml", "w", encoding="utf-8", newline="\\n") as handle:
        handle.write('"$task.process":\\n')
        handle.write("    gwaslab: {}\\n".format(gl.__version__))
        handle.write("    python: {}\\n".format(sys.version.split()[0]))


if __name__ == "__main__":
    sys.exit(main())
