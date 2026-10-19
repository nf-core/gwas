#!/usr/bin/env python3
"""Prepare one analysis unit's phenotype and covariate files for downstream programmes.

The three association and heritability programmes this pipeline drives disagree irreconcilably about
two things. REGENIE requires a header row and GCTA forbids one, so the phenotype and the covariates
are written twice, once with a header and once without. And each takes a different trait selector --
a column name, or a one-based trait index -- so the trait is always written to the third column of a
three-column file under the constant name PHENO, which turns GCTA's and LDAK's `--mpheno` into the
constant 1 and makes REGENIE's phenotype-named output filename deterministic.

At most eight files are written, all tab-delimited with LF line endings:

    <prefix>.pheno              FID IID PHENO                          always
    <prefix>.qcovar             FID IID <quantitative names>           when the row supplied them
    <prefix>.catcovar           FID IID <categorical names>            when the row supplied them
    <prefix>.covar              FID IID <quantitative> <categorical>   when the row supplied either
    <prefix>.noheader.pheno     the same content minus line 1          always
    <prefix>.noheader.qcovar    the same content minus line 1          when the row supplied them
    <prefix>.noheader.catcovar  the same content minus line 1          when the row supplied them
    <prefix>.adjustcovar        quantitative plus encoded factors       when the row supplied either

Standard library only, deliberately: the inputs run to a few thousand rows at most, so a dataframe
dependency would add container weight and a second version to report for no gain.
"""

import sys

# Interpolated as JSON string literals, which are valid Python syntax. Optional paths are
# serialized as empty strings when the row did not supply them.
PHENOTYPE_FILE = ${phenotype_literal}
QUANT_COVARIATES_FILE = ${quant_covariates_literal}
CAT_COVARIATES_FILE = ${cat_covariates_literal}
PHENOTYPE_COLUMN = ${phenotype_column_literal}
TRAIT_TYPE = ${trait_type_literal}
CASE_VALUE = ${case_value_literal}
CONTROL_VALUE = ${control_value_literal}
COVARIATE_COMPLETENESS = ${covariate_completeness_literal}
PREFIX = ${prefix_literal}
ANALYSIS_ID = ${analysis_id_literal}
PROCESS_NAME = ${task_process_literal}

MISSING = "NA"

# Keep raw-value errors useful without allowing an unbounded message for bad input.
MAX_UNMATCHED_VALUES = 10

# `NA` is the one missing code all three programmes accept. `-9` is read as missing on the way in
# because PLINK and GCTA both write it, but it is never written out: LDAK and REGENIE reject it, and
# PLINK 2 (which the pipeline no longer feeds directly) errors out when a `-9` shares a file with a
# value in (-9, 10]. The empty string is a genuine member: it is how a tab-delimited file spells a
# missing cell, and `split_row` preserves it.
MISSING_TOKENS = frozenset(["", "na", "nan", "-9"])

def fail(message):
    """Abort naming the analysis unit, so a large samplesheet can be fixed without bisecting it."""
    sys.exit("[nf-core/gwas] ERROR: analysis '{}': {}".format(ANALYSIS_ID, message))


def split_row(line):
    """Split one input line into fields, per line, on whichever delimiter it actually uses.

    PLINK writes and accepts both tab- and space-delimited files and a researcher's own file may be
    either, so neither delimiter can simply be assumed. Splitting on arbitrary whitespace alone is
    wrong, though: it collapses a run of tabs, so an empty cell -- the ordinary way a tab-delimited
    file spells a missing value -- would vanish and leave the row looking ragged rather than missing.
    A line carrying a tab is therefore a tab-delimited line and is split on tabs, which preserves the
    empty cell for `is_missing` to normalise to NA. No value is ever split by either branch: the
    samplesheet contract forbids whitespace inside an identifier or a trait value.
    """
    line = line.rstrip("\\r\\n")
    return line.split("\\t") if "\\t" in line else line.split()


def read_table(path, role):
    """Read one FID/IID-keyed input file, rejecting the shapes that cannot be normalised."""
    with open(path) as handle:
        rows = [split_row(line) for line in handle]
    # A row of nothing but empty cells is a blank line, not a sample with a missing everything.
    rows = [row for row in rows if any(field.strip() for field in row)]
    if not rows:
        fail("the {} file '{}' is empty".format(role, path))

    header = rows[0]
    body = rows[1:]
    if [name.upper() for name in header[:2]] != ["FID", "IID"]:
        fail(
            "the {} file '{}' must start with an FID column and an IID column, but starts with {}".format(
                role, path, ", ".join("'{}'".format(name) for name in header[:2]) or "nothing"
            )
        )

    for number, row in enumerate(body, start=2):
        if len(row) != len(header):
            fail(
                "the {} file '{}' declares {} columns but line {} carries {} fields".format(
                    role, path, len(header), number, len(row)
                )
            )

    seen = set()
    for row in body:
        identity = (row[0], row[1])
        if identity in seen:
            fail("the {} file '{}' declares sample '{} {}' more than once".format(role, path, row[0], row[1]))
        seen.add(identity)

    return header, body


def is_missing(value):
    return value.strip().lower() in MISSING_TOKENS


def normalise_trait(value):
    """Recode one source trait value to the coding every downstream programme accepts.

    Binary traits become 0/1/NA -- the only coding REGENIE, GCTA and LDAK all read the same way --
    by comparing the source cell against the declared case and control values as strings. The
    samplesheet carries both as text precisely so that PLINK's 1/2, a 0/1 file and labels such as
    'Case' all work without the pipeline guessing which convention is in force. A cell matching
    neither is missing. The ingress checks below require both declared binary values to occur.
    """
    if TRAIT_TYPE == "binary":
        stripped = value.strip()
        if stripped == CASE_VALUE:
            return "1"
        if stripped == CONTROL_VALUE:
            return "0"
        return MISSING
    if is_missing(value):
        return MISSING
    try:
        float(value)
    except ValueError:
        return MISSING
    return value


def normalise_covariate_row(row):
    """Covariates pass through verbatim; only their missing code is rewritten."""
    return row[:2] + [MISSING if is_missing(value) else value for value in row[2:]]


def load_covariates(path, role):
    """Read and normalise one covariate file, or return None when the row supplied none."""
    if not path:
        return None
    header, body = read_table(path, role)
    return header, [normalise_covariate_row(row) for row in body]


def merge_covariates(quant, cat):
    """Build the single file REGENIE takes, quantitative columns first.

    The join is an inner one on FID and IID: a sample described by only one of the two files has an
    incomplete covariate vector and would be dropped by every consumer anyway.
    """
    if quant is None:
        return cat
    if cat is None:
        return quant
    quant_header, quant_body = quant
    cat_header, cat_body = cat
    cat_by_identity = {(row[0], row[1]): row[2:] for row in cat_body}
    body = [row + cat_by_identity[(row[0], row[1])] for row in quant_body if (row[0], row[1]) in cat_by_identity]
    return quant_header + cat_header[2:], body


# --adjust-grm rejects --factors on 6.1, 6.2 and
# 6.3, so the categorical adjustment design is treatment-coded here. Retire when an LDAK release accepts
# --factors with --adjust-grm.
def adjustment_covariates(quant, cat):
    """Build the numerical design accepted by LDAK ``--adjust-grm``.

    LDAK's estimators accept categorical covariates through ``--factors``, but ``--adjust-grm``
    rejects that flag. The equivalent matrix-adjustment design therefore carries quantitative
    columns unchanged and treatment-codes every categorical column against its lexically first
    observed level. Missing factor values become missing across every dummy for that factor.
    """
    if cat is None:
        return quant

    cat_header, cat_body = cat
    factor_levels = []
    adjustment_header = ["FID", "IID"]
    for index, name in enumerate(cat_header[2:], start=2):
        levels = sorted({row[index] for row in cat_body if not is_missing(row[index])})
        encoded_levels = levels[1:]
        factor_levels.append((index, encoded_levels))
        adjustment_header.extend("{}_{}".format(name, level) for level in encoded_levels)

    encoded_cat_by_identity = {}
    for row in cat_body:
        encoded = []
        for index, levels in factor_levels:
            value = row[index]
            encoded.extend(
                [MISSING] * len(levels)
                if is_missing(value)
                else ["1" if value == level else "0" for level in levels]
            )
        encoded_cat_by_identity[(row[0], row[1])] = encoded

    if quant is None:
        body = [row[:2] + encoded_cat_by_identity[(row[0], row[1])] for row in cat_body]
        return adjustment_header, body

    quant_header, quant_body = quant
    body = [
        row + encoded_cat_by_identity[(row[0], row[1])]
        for row in quant_body
        if (row[0], row[1]) in encoded_cat_by_identity
    ]
    return quant_header + adjustment_header[2:], body


def write_lines(path, lines):
    with open(path, "w", newline="\\n") as handle:
        handle.writelines(line + "\\n" for line in lines)
        handle.flush()


def write_table(suffix, header, body):
    """Write the headered serialisation REGENIE takes and the headerless one GCTA and LDAK take."""
    rows = ["\\t".join(row) for row in body]
    write_lines("{}.{}".format(PREFIX, suffix), ["\\t".join(header)] + rows)
    write_lines("{}.noheader.{}".format(PREFIX, suffix), rows)


def raw_value_counts(raw_values):
    """Count source-cell spellings before any normalisation can hide them."""
    counts = {}
    for value in raw_values:
        counts[value] = counts.get(value, 0) + 1
    return counts


def format_raw_value_examples(counts, kind):
    """Render bounded, deterministic source-value examples without an unbounded diagnostic."""
    if not counts:
        return "none"
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    displayed = ranked[:MAX_UNMATCHED_VALUES]
    summary = ", ".join("'{}' ({} samples)".format(value, count) for value, count in displayed)
    elided = len(ranked) - len(displayed)
    if elided:
        summary += "; {} further distinct {} values omitted".format(elided, kind)
    return summary


def binary_match_counts(raw_values):
    """Count declared binary values with the exact stripping semantics normalise_trait uses."""
    case_matches = 0
    control_matches = 0
    for value in raw_values:
        stripped = value.strip()
        if stripped == CASE_VALUE:
            case_matches += 1
        if stripped == CONTROL_VALUE:
            control_matches += 1
    return case_matches, control_matches


# LDAK 6.1 (pinned genomedk build) fits a stale read-buffer
# value for a missing --covar cell at exit 0. Retire when the LDAK analysis-row modules move to 6.3,
# which mean-imputes as documented.
def missing_covariate_cells(covariates, phenotyped):
    """List every (FID, IID, column) whose covariate cell is missing for a sample that has a trait value.

    Samples the estimator would drop for a missing phenotype are not policed: their covariates are never read.
    """
    if covariates is None:
        return []
    header, body = covariates
    found = []
    for row in body:
        if (row[0], row[1]) not in phenotyped:
            continue
        for index, name in enumerate(header[2:], start=2):
            if is_missing(row[index]):
                found.append((row[0], row[1], name))
    return found


def numeric_match_count(raw_values):
    """Count source cells that pass the same numeric parsing normalise_trait applies."""
    matches = 0
    for value in raw_values:
        if is_missing(value):
            continue
        try:
            float(value)
        except ValueError:
            continue
        matches += 1
    return matches


phenotype_header, phenotype_body = read_table(PHENOTYPE_FILE, "phenotype")
if PHENOTYPE_COLUMN not in phenotype_header:
    fail(
        "phenotype column '{}' is not present in '{}', which declares {}".format(
            PHENOTYPE_COLUMN, PHENOTYPE_FILE, ", ".join("'{}'".format(name) for name in phenotype_header)
        )
    )
trait_index = phenotype_header.index(PHENOTYPE_COLUMN)
if trait_index < 2:
    fail(
        "phenotype column '{}' is the {} identifier column of '{}', not a trait".format(
            PHENOTYPE_COLUMN, "family" if trait_index == 0 else "individual", PHENOTYPE_FILE
        )
    )

raw_trait_values = [row[trait_index] for row in phenotype_body]
raw_counts = raw_value_counts(raw_trait_values)
case_raw_matches, control_raw_matches = binary_match_counts(raw_trait_values)
numeric_raw_matches = numeric_match_count(raw_trait_values)
normalised_trait_values = [normalise_trait(value) for value in raw_trait_values]
trait_rows = [
    [row[0], row[1], normalised_value]
    for row, normalised_value in zip(phenotype_body, normalised_trait_values)
]

quant_covariates = load_covariates(QUANT_COVARIATES_FILE, "quantitative covariate")
cat_covariates = load_covariates(CAT_COVARIATES_FILE, "categorical covariate")

# The merged file has no headerless counterpart: GCTA and LDAK take the two kinds of covariate
# through separate flags and never see a merged file.
merged = merge_covariates(quant_covariates, cat_covariates)
adjustment = adjustment_covariates(quant_covariates, cat_covariates)

if TRAIT_TYPE == "binary":
    if case_raw_matches == 0 or control_raw_matches == 0:
        fail(
            "phenotype source '{}' column '{}' declares case_value '{}' and control_value '{}'; "
            "case_value '{}' matched {}; control_value '{}' matched {}; {}".format(
                PHENOTYPE_FILE,
                PHENOTYPE_COLUMN,
                CASE_VALUE,
                CONTROL_VALUE,
                CASE_VALUE,
                case_raw_matches,
                CONTROL_VALUE,
                control_raw_matches,
                "raw values: {}".format(format_raw_value_examples(raw_counts, "raw")),
            )
        )
elif TRAIT_TYPE == "quantitative" and numeric_raw_matches == 0:
    fail(
        "phenotype source '{}' column '{}' matched 0 numeric values; raw values: {}".format(
            PHENOTYPE_FILE, PHENOTYPE_COLUMN, format_raw_value_examples(raw_counts, "raw")
        )
    )

# LDAK 6.1 (pinned genomedk build) fits a stale read-buffer
# value for a missing --covar cell at exit 0. Retire when the LDAK analysis-row modules move to 6.3,
# which mean-imputes as documented.
if COVARIATE_COMPLETENESS == "required":
    phenotyped = set((row[0], row[1]) for row in trait_rows if row[2] != MISSING)
    incomplete = missing_covariate_cells(quant_covariates, phenotyped) + missing_covariate_cells(
        cat_covariates, phenotyped
    )
    if incomplete:
        affected = set((fid, iid) for fid, iid, _name in incomplete)
        displayed = sorted(incomplete)[:MAX_UNMATCHED_VALUES]
        elided = len(incomplete) - len(displayed)
        fail(
            "selects method(s) that read covariates through LDAK --covar/--factors, which treat a missing "
            "cell as a value rather than excluding the sample; {} sample(s) have missing covariate cells "
            "(first {}: {}{}). Remove those samples from the phenotype file or complete the covariates".format(
                len(affected),
                len(displayed),
                ", ".join("{} {} {}".format(fid, iid, name) for fid, iid, name in displayed),
                "; {} further missing cells omitted".format(elided) if elided else "",
            )
        )

write_table("pheno", ["FID", "IID", "PHENO"], trait_rows)
if quant_covariates is not None:
    write_table("qcovar", quant_covariates[0], quant_covariates[1])
if cat_covariates is not None:
    write_table("catcovar", cat_covariates[0], cat_covariates[1])
if merged is not None:
    write_lines("{}.covar".format(PREFIX), ["\\t".join(row) for row in [merged[0]] + merged[1]])
if adjustment is not None:
    write_lines("{}.adjustcovar".format(PREFIX), ["\\t".join(row) for row in adjustment[1]])

# Written here rather than captured by an `eval` output, which Nextflow allows only on a Bash script.
write_lines("versions.yml", ['"{}":'.format(PROCESS_NAME), "    python: {}".format(sys.version.split()[0])])
