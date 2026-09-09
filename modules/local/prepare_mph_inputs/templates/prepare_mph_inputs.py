#!/usr/bin/env python3
"""Write exactly the CSVs MPH consumes for one analysis unit.

MPH's file interface differs from every other individual-level estimator this pipeline drives, in four ways:

    it keys the phenotype and covariate tables by IID alone, while indexing the relationship
        matrix itself by the *order* of that matrix's `.grm.iid`. The adapter follows those native
        identities and that declared order. It cannot prove that a `.grm.iid` still matches its own
        `.grm.bin`; the route avoids that residual hazard by building both from the same genotype view.
    an empty field is its only missing representation, and every other cell is parsed as a number. This is
        not a preference: the pipeline's three missing spellings are therefore rewritten to empty fields.
    it synthesises an intercept only when no covariate is named. A covariate file given with
        `--covariate_names` fits exactly the named columns, so a covariate-adjusted fit without an explicit
        column of ones is a no-intercept model. This adapter always writes that column and names it first.
    it never expands a categorical covariate. A factor column passed through verbatim becomes a
        numeric covariate whose levels are read as magnitudes, and the resulting rank deficiency is only a
        warning at exit 0. Categorical columns are therefore dummy-encoded here, by the same rule
        `prepare_phenotype_inputs` applies for LDAK's matrix-adjustment design, so the two encodings agree.

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

# The pipeline's own missing spellings, matching `prepare_phenotype_inputs`. Every one of them becomes an empty
# CSV field, which is MPH's only missing representation.
MISSING_TOKENS = frozenset(["", "na", "nan", "-9"])


def fail(message):
    sys.exit("[nf-core/gwas] ERROR: analysis '{}': {}".format(ANALYSIS_ID, message))


def is_missing(value):
    return value.strip().lower() in MISSING_TOKENS


def split_row(line):
    """Split one input line on whichever delimiter it actually uses, per line.

    The rule is `prepare_phenotype_inputs`' and it is copied rather than approximated, because getting it
    wrong here is a silent row shift rather than a parse error. Splitting on arbitrary whitespace collapses a
    run of tabs, so a genuinely empty cell -- the ordinary way a tab-delimited file spells a missing value --
    would vanish and the values to its right would each move one column left, arriving at MPH under the wrong
    covariate name. A line carrying a tab is therefore a tab-delimited line and is split on tabs.
    """
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
    """Read a prepared covariate table: two identifier columns, then one column per named covariate.

    The shape of the identifier columns is not re-checked here. Ingress already admitted this file and
    `prepare_phenotype_inputs` wrote it, so a check would only be able to fail a row the pipeline had already
    accepted -- and the incoming header is never forwarded to MPH in any case, because this module emits its
    own `IID,...` line from the covariate names it derives.
    """
    rows = read_table(path, role)
    header = rows[0][1]
    if len(header) < 3:
        fail(
            "the {} file '{}' has {} columns; two identifier columns and at least one named covariate "
            "column are required".format(role, path, len(header))
        )
    # A ragged row would otherwise be handed to MPH as a short CSV line under a full-width header, which MPH
    # reads without complaint.
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

    # MPH indexes the matrix by this file's line order and never cross-checks it: a `.grm.iid` reordered
    # against its own `.grm.bin` gives a complete, stable, wrong estimate at exit 0 (pve 0.1156 -> -0.2995,
    # 5/5 runs, pinned image). Not asserted here: MPH_MAKEGRM builds every matrix from this view.
    return fam_by_iid, grm_order


def read_traits(fam_by_iid):
    """Read the headerless `FID IID trait...` table the pipeline prepares, keyed by IID as MPH keys it."""
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
    """Reapply the pipeline's covariate encoding rule against the headered prepared tables.

    MPH's covariate interface is name-keyed, so this adapter reads the headered `.qcovar`/`.catcovar` tables
    rather than the headerless numeric design LDAK consumes. The encoding itself is the rule of
    `prepare_phenotype_inputs.adjustment_covariates`: quantitative columns verbatim, every categorical column
    treatment-coded against its lexically first observed level, and a missing factor value missing across all
    of that factor's dummies. Keeping the two in step is what a parity test exists for.
    """
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
    # The inner join `prepare_phenotype_inputs` applies for the same reason: a sample described by only one of
    # the two files has an incomplete covariate vector, which MPH would drop anyway.
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
