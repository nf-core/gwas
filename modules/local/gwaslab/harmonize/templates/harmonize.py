#!/usr/bin/env python3

import json
import math
import numbers
import os
import shlex
import sys

# The packaged image is read-only at the locations selected by Numba and Matplotlib.
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.abspath(".numba_cache"))
os.environ.setdefault("MPLCONFIGDIR", os.path.abspath(".matplotlib"))

import gwaslab as gl


def optional_path(value: str):
    return None if value in {"", "[]", "null"} else value


def load_format(value: str):
    """Return a GWASLab format name or explicit constructor column mapping."""
    if not value.lstrip().startswith("{"):
        if value.startswith("auto"):
            raise ValueError("Automatic input-format detection is not supported")
        return value, {}
    mapping = json.loads(value)
    if not isinstance(mapping, dict) or not mapping:
        raise ValueError("An explicit input column mapping must be a non-empty JSON object")
    return None, mapping


# GWASLab writes every effect size, its confidence bounds and every test statistic with a fixed
# four-decimal format, which is lossy at the magnitudes a well-powered GWAS produces: a common
# variant at N=450,000 with BETA=3.12e-05 and SE=9.8e-06 is written as `0.0000  0.0000`, destroying
# the effect and making an inverse-variance weight `1/SE**2` divide by zero downstream. Six-digit
# scientific notation keeps the same significant digits at every magnitude, and comfortably exceeds
# the six significant digits the association programmes themselves write. GWASLab's remaining
# defaults are left alone: allele frequencies already use a significant-digit format and P already
# uses scientific notation.
FLOAT_FORMATS = {
    column: "{:.6e}"
    for column in [
        "BETA",
        "SE",
        "BETA_95L",
        "BETA_95U",
        "OR",
        "OR_95L",
        "OR_95U",
        "HR",
        "HR_95L",
        "HR_95U",
        "Z",
        "CHISQ",
        "F",
        "MLOG10P",
    ]
}
prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
args = shlex.split("$task.ext.args" if "$task.ext.args" != "null" else "")
allowed = {"--keep-invalid", "--ref-alt-freq"}
unknown = [arg for arg in args if arg.startswith("--") and arg not in allowed]
if unknown:
    raise ValueError(f"Unsupported option(s): {', '.join(unknown)}")

remove_invalid = "--keep-invalid" not in args
ref_alt_freq = None
if "--ref-alt-freq" in args:
    index = args.index("--ref-alt-freq")
    if index + 1 == len(args):
        raise ValueError("--ref-alt-freq requires a value")
    ref_alt_freq = args[index + 1]

builds = {"GRCh37": "19", "GRCh38": "38"}
if "$genome_build" not in builds:
    raise ValueError("genome_build must be GRCh37 or GRCh38")

fmt, mapping = load_format(r'''$input_format''')
sumstats = gl.Sumstats(
    "$sumstats",
    fmt=fmt,
    build=builds["$genome_build"],
    species="homo sapiens",
    **mapping,
)
sumstats.harmonize(
    ref_seq=optional_path("$reference_fasta"),
    ref_rsid_vcf=optional_path("$rsid_reference_vcf"),
    ref_infer=optional_path("$strand_reference_vcf"),
    ref_alt_freq=ref_alt_freq,
    threads=$task.cpus,
    remove=remove_invalid,
    fix_id_kwargs={"fixchrpos": True},
    sweep_mode=False,
)
if "$meta.method" == "ldak_kvik":
    if "N_EFF" not in sumstats.data.columns:
        raise ValueError("LDAK-KVIK harmonisation requires GWASLab's N_EFF column")
    effective_n = sumstats.data["N_EFF"]
    invalid_effective_n = effective_n.isna() | ~effective_n.map(
        lambda value: isinstance(value, numbers.Real) and math.isfinite(value) and value > 0
    )
    if invalid_effective_n.any():
        raise ValueError(
            "LDAK-KVIK effective analysis sizes must be numeric, finite, non-null and greater than zero"
        )
    if "N" in sumstats.data.columns:
        raise ValueError("LDAK-KVIK harmonisation cannot promote N_EFF because N already exists")
    sumstats.data.rename(columns={"N_EFF": "N"}, inplace=True)
sumstats.to_format(
    prefix,
    fmt="gwaslab",
    tab_fmt="tsv",
    gzip=True,
    output_log=True,
    float_formats=FLOAT_FORMATS,
)

with open("versions.yml", "w", encoding="utf-8") as versions:
    versions.write(f'"$task.process":\\n')
    versions.write(f"    gwaslab: {gl.__version__}\\n")
    versions.write(f"    python: {sys.version.split()[0]}\\n")
