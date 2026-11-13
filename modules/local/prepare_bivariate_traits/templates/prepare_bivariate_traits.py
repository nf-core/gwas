#!/usr/bin/env python3
"""Build the ordered two-trait GCTA phenotype and relationship-owned covariates."""

import sys

LEFT_PHENOTYPE = ${left_phenotype_literal}
RIGHT_PHENOTYPE = ${right_phenotype_literal}
PAIR_QUANT_COVARIATES = ${pair_quant_covariates_literal}
PAIR_CAT_COVARIATES = ${pair_cat_covariates_literal}
PREFIX = ${prefix_literal}
REQUEST_ID = ${request_id_literal}
PROCESS_NAME = ${task_process_literal}

MISSING = "NA"
MISSING_TOKENS = frozenset(["", "na", "nan", "-9"])


def fail(message):
    sys.exit("[nf-core/gwas] ERROR: pair request '{}': {}".format(REQUEST_ID, message))


def split_row(line):
    line = line.rstrip("\\r\\n")
    return line.split("\\t") if "\\t" in line else line.split()


def is_missing(value):
    return value.strip().lower() in MISSING_TOKENS


def read_phenotype(path):
    values = {}
    with open(path) as handle:
        for line in handle:
            row = split_row(line)
            values[(row[0], row[1])] = row[2]
    return values


def read_covariates(path, role):
    if not path:
        return None, []
    with open(path) as handle:
        rows = [split_row(line) for line in handle]
    rows = [row for row in rows if row and any(field.strip() for field in row)]
    if not rows:
        fail("the {} file '{}' is empty".format(role, path))
    header = rows[0]
    if len(header) < 3 or [value.upper() for value in header[:2]] != ["FID", "IID"]:
        fail("the {} file '{}' must begin with FID, IID and at least one covariate".format(role, path))
    if len(set(header)) != len(header):
        fail("the {} file '{}' contains duplicate column names".format(role, path))
    body = []
    seen = set()
    for line_number, row in enumerate(rows[1:], start=2):
        if len(row) != len(header):
            fail("the {} file '{}' declares {} columns but line {} contains {} fields".format(role, path, len(header), line_number, len(row)))
        identity = (row[0], row[1])
        if identity in seen:
            fail("the {} file '{}' declares sample '{} {}' more than once".format(role, path, *identity))
        seen.add(identity)
        body.append(row[:2] + [MISSING if is_missing(value) else value for value in row[2:]])
    if not body:
        fail("the {} file '{}' contains no sample rows".format(role, path))
    return header, body


def write_rows(path, rows):
    with open(path, "w", newline="") as handle:
        for row in rows:
            handle.write("\\t".join(row) + "\\n")


left = read_phenotype(LEFT_PHENOTYPE)
right = read_phenotype(RIGHT_PHENOTYPE)
union = sorted(set(left) | set(right))
phenotype_rows = [
    [fid, iid, left.get((fid, iid), MISSING), right.get((fid, iid), MISSING)]
    for fid, iid in union
]
write_rows(PREFIX + ".pheno", phenotype_rows)

quant_header, quant_rows = read_covariates(PAIR_QUANT_COVARIATES, "relationship quantitative covariate")
cat_header, cat_rows = read_covariates(PAIR_CAT_COVARIATES, "relationship categorical covariate")
# The headered files carry the header this module has already validated, so the two serialisations are the same
# rows and differ only in whether the column names survive. GCTA rejects a header; MPH addresses its covariates
# by name.
if quant_header:
    write_rows(PREFIX + ".qcovar", [quant_header] + quant_rows)
    write_rows(PREFIX + ".noheader.qcovar", quant_rows)
if cat_header:
    write_rows(PREFIX + ".catcovar", [cat_header] + cat_rows)
    write_rows(PREFIX + ".noheader.catcovar", cat_rows)

with open("versions.yml", "w", newline="") as handle:
    handle.write('"{}":\\n'.format(PROCESS_NAME))
    handle.write("    python: {}\\n".format(sys.version.split()[0]))
