#!/usr/bin/env python3
"""Serialize FID/IID trait and covariate tables as IID-keyed MPH CSVs.

Trait columns retain the caller's explicit order. Quantitative columns pass through;
factors are treatment-coded against their lexically first observed level. Covariate
files receive an explicit intercept. The caller declares missing-value spellings.
The supplied sample map resolves FIDs, and the native matrix IID list orders rows.
"""

import json
import sys

PHENOTYPE_TABLE = $phenotype_table_literal
QUANT_COVARIATES = $quant_covariates_literal
CAT_COVARIATES = $cat_covariates_literal
GRM_IID = $grm_iid_literal
FAM = $fam_literal
PREFIX = $prefix_literal
ANALYSIS_ID = $analysis_id_literal
PROCESS_NAME = $task_process_literal
TRAIT_NAMES = json.loads($trait_names_literal)

# Missing-value spellings are supplied by the caller, never inferred.
MISSING_TOKENS = frozenset(value.strip().lower() for value in json.loads($missing_tokens_literal))


def fail(message):
    sys.exit("prepare_mph_inputs: sample set '{}': {}".format(ANALYSIS_ID, message))


def is_missing(value):
    return value.strip().lower() in MISSING_TOKENS


def split_row(line):
    """Preserve empty tab-delimited cells; accept whitespace-delimited nonempty fields."""
    line = line.rstrip("\\r\\n")
    return line.split("\\t") if "\\t" in line else line.split()


def read_table(path, role):
    rows = []
    with open(path) as handle:
        for number, line in enumerate(handle, start=1):
            fields = split_row(line)
            # A row of nothing but empty cells is a blank line, not a sample missing everything.
            if any(field.strip() for field in fields):
                rows.append((number, fields))
    if not rows:
        fail("the {} file '{}' has no data rows".format(role, path))
    return rows


def read_headered_table(path, role):
    """Read a headered FID IID covariate table."""
    rows = read_table(path, role)
    header = rows[0][1]
    if len(header) < 3:
        fail(
            "the {} file '{}' has {} columns; two identifier columns and at least one named covariate "
            "column are required".format(role, path, len(header))
        )
    # A ragged row cannot be assigned to columns -- which field belongs to which name is undecidable -- so it
    # is named and refused here rather than serialised.
    for number, fields in rows[1:]:
        if len(fields) != len(header):
            fail(
                "the {} file '{}' declares {} columns but line {} carries {} fields".format(
                    role, path, len(header), number, len(fields)
                )
            )
    return header, [fields for _number, fields in rows[1:]]


def write_lines(path, lines):
    with open(path, "w", newline="\\n") as handle:
        handle.writelines(line + "\\n" for line in lines)
        handle.flush()


def read_sample_order():
    """Read the genotype identities and the matrix order that MPH applies to the prepared tables."""
    fam_rows = read_table(FAM, "genotype FAM")
    fam_by_iid = {}
    for number, fields in fam_rows:
        if len(fields) < 2:
            fail("row {} of the genotype FAM '{}' has fewer than two columns".format(number, FAM))
        fid, iid = fields[0], fields[1]
        fam_by_iid[iid] = fid

    with open(GRM_IID) as handle:
        grm_order = [line.rstrip("\\r\\n") for line in handle if line.rstrip("\\r\\n") != ""]

    return fam_by_iid, grm_order


def read_traits(fam_by_iid):
    """Read the headerless FID IID trait table, keyed by IID."""
    expected_columns = 2 + len(TRAIT_NAMES)
    rows = read_table(PHENOTYPE_TABLE, "phenotype")
    traits_by_iid = {}
    for number, fields in rows:
        if len(fields) != expected_columns:
            fail(
                "row {} of the phenotype table '{}' has {} columns, expected {} for traits {}".format(
                    number, PHENOTYPE_TABLE, len(fields), expected_columns, ", ".join(TRAIT_NAMES)
                )
            )
        fid, iid = fields[0], fields[1]
        if iid not in fam_by_iid:
            continue
        if iid in traits_by_iid:
            fail("the phenotype table '{}' lists IID '{}' more than once".format(PHENOTYPE_TABLE, iid))
        traits_by_iid[iid] = ["" if is_missing(value) else value for value in fields[2:]]
    return traits_by_iid


def encode_covariates():
    """Treatment-code categorical columns and join quantitative/categorical designs on FID/IID."""
    quant = read_headered_table(QUANT_COVARIATES, "quantitative covariate") if QUANT_COVARIATES else None
    cat = read_headered_table(CAT_COVARIATES, "categorical covariate") if CAT_COVARIATES else None
    if quant is None and cat is None:
        return None

    quant_names = []
    quant_by_identity = {}
    if quant is not None:
        quant_header, quant_body = quant
        quant_names = quant_header[2:]
        for row in quant_body:
            values = ["" if is_missing(value) else value for value in row[2:]]
            quant_by_identity[(row[0], row[1])] = values

    dummy_names = []
    cat_by_identity = {}
    if cat is not None:
        cat_header, cat_body = cat
        factor_levels = []
        for index, name in enumerate(cat_header[2:], start=2):
            levels = sorted({row[index] for row in cat_body if not is_missing(row[index])})
            encoded_levels = levels[1:]
            factor_levels.append((index, encoded_levels))
            dummy_names.extend("{}_{}".format(name, level) for level in encoded_levels)
        for row in cat_body:
            values = []
            for index, levels in factor_levels:
                value = row[index]
                if is_missing(value):
                    values.extend([""] * len(levels))
                else:
                    values.extend("1" if value == level else "0" for level in levels)
            cat_by_identity[(row[0], row[1])] = values

    if cat is None:
        return quant_names, quant_by_identity
    if quant is None:
        return dummy_names, cat_by_identity
    # Only identities represented by both sources have a complete joined design.
    merged = {
        identity: quant_by_identity[identity] + cat_by_identity[identity]
        for identity in quant_by_identity
        if identity in cat_by_identity
    }
    return quant_names + dummy_names, merged


def main():
    fam_by_iid, grm_order = read_sample_order()
    traits_by_iid = read_traits(fam_by_iid)
    covariates = encode_covariates()

    covariate_names = []
    covariate_by_identity = {}
    if covariates is not None:
        source_names, covariate_by_identity = covariates
        # `intercept` is written first and is always 1. MPH synthesises an intercept only when no covariate is
        # named, so omitting this column would silently fit a model through the origin.
        covariate_names = ["intercept"] + list(source_names)

    written = [iid for iid in grm_order if iid in traits_by_iid]

    phenotype_lines = ["IID," + ",".join(TRAIT_NAMES)]
    covariate_lines = ["IID," + ",".join(covariate_names)] if covariates is not None else None

    for iid in written:
        values = traits_by_iid[iid]
        phenotype_lines.append(",".join([iid] + values))

        if covariates is not None:
            row = covariate_by_identity.get((fam_by_iid[iid], iid))
            # A sample the covariate tables do not describe is written with empty cells rather than omitted, so
            # the two files stay row-aligned.
            row = [""] * (len(covariate_names) - 1) if row is None else list(row)
            covariate_lines.append(",".join([iid, "1"] + row))

    write_lines("{}.mph.pheno.csv".format(PREFIX), phenotype_lines)
    if covariate_lines is not None:
        write_lines("{}.mph.covar.csv".format(PREFIX), covariate_lines)

    # Written here rather than captured by an `eval` output, which Nextflow allows only on a Bash script.
    write_lines("versions.yml", ['"{}":'.format(PROCESS_NAME), "    python: {}".format(sys.version.split()[0])])
    return 0


if __name__ == "__main__":
    sys.exit(main())
