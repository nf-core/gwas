#!/usr/bin/env python3

import json
import os
import sys

# Importing GWASLab 4.1.9 as a non-root container user raises RuntimeError when
# Numba cache=True targets the read-only package directory. Retire when the bioconda::gwaslab pin moves past the fix.
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.abspath(".numba_cache"))
os.environ.setdefault("MPLCONFIGDIR", os.path.abspath(".matplotlib"))

import gwaslab as gl


def optional_path(value: str):
    return None if value in {"", "[]", "null"} else value


def load_format(value: str):
    """Return a GWASLab format name or explicit constructor column mapping."""
    if not value.lstrip().startswith("{"):
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
# GWASLab 4.1.9 to_format
# defaults to {:.4f}, which zeroes small BETA/SE.
# Retire if float_formats becomes GWASLab's default.
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
    threads=$task.cpus,
)
if "$meta.method" == "ldak_kvik":
    if "N_EFF" not in sumstats.data.columns:
        raise ValueError("LDAK-KVIK harmonisation requires GWASLab's N_EFF column")
    if "N" in sumstats.data.columns:
        raise ValueError("LDAK-KVIK harmonisation cannot promote N_EFF because N already exists")
    sumstats.data.rename(columns={"N_EFF": "N"}, inplace=True)
if "P" not in sumstats.data.columns:
    # REGENIE and some supported external formats report -log10(P) instead of P. The
    # pipeline-standard representation requires P, so ask GWASLab to derive it from whichever
    # standardised test-statistic column the explicit input format supplied.
    sumstats.fill_data(to_fill=["P"])
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
