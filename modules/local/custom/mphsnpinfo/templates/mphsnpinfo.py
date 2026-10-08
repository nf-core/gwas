#!/usr/bin/env python3
"""Convert ordered disjoint named SNP groups to native MPH weight columns."""

import csv
import sys

bim_path = $bim_literal
prefix = $prefix_literal
manifest_path = $manifest_literal
columns = $columns_literal
process_name = $task_process_literal


def fail(message):
    sys.exit("MPH SNP information: {}".format(message))


# Validates the caller's column names: `--snp_weight_name` selects a column by name, so a repeated name
# would leave MPH silently building the matrix from whichever duplicate it matched first.
if len(set(columns)) != len(columns):
    fail("weight column names repeat: {}".format(columns))

with open(bim_path) as handle:
    variants = [line.split()[1] for line in handle]

known_variants = set(variants)
membership = {}

# Input conversion is required because MPH consumes numeric weight columns rather than SNP lists.
with open(manifest_path) as handle:
    records = list(csv.DictReader(handle, delimiter="\\t"))
for record in records:
    # Membership is keyed by the column name the row declares, not by the row's position, so the written
    # columns follow the caller's weight column order whatever order the manifest rows arrive in.
    index = columns.index(record["group_id"])
    with open(record["group_filename"]) as handle:
        group = handle.read().split()
    for name in group:
        # Validates a user-supplied resource file: a group list may only name variants of the supplied BIM.
        if name not in known_variants:
            fail("SNP '{}' is in group '{}' but not in the BIM".format(name, record["group_id"]))
        # Validates a user-supplied resource file: overlapping groups would give one variant weight in two
        # components, which is not a partition and not what the written CSV could represent.
        if name in membership:
            fail("SNP '{}' is in more than one group file; groups must be disjoint".format(name))
        membership[name] = index

counts = [0] * len(columns)
unassigned = 0
with open("{}.snp_info.csv".format(prefix), "w", newline="\\n") as handle:
    writer = csv.writer(handle, lineterminator="\\n")
    writer.writerow(["SNP"] + columns)
    for name in variants:
        index = membership.get(name)
        if index is None:
            unassigned += 1
        else:
            counts[index] += 1
        writer.writerow([name] + ["1" if index == column else "0" for column in range(len(columns))])

# `__unassigned__` is the number of BIM rows MPH is told to ignore and `__bim_rows__` the number it will
# report as found, so the pair reconciles the supplied groups against the BIM rows.
with open("{}.snp_info.counts.tsv".format(prefix), "w", newline="\\n") as handle:
    handle.write("group_id\\tpredictor_count\\n")
    for name, count in zip(columns, counts):
        handle.write("{}\\t{}\\n".format(name, count))
    handle.write("__unassigned__\\t{}\\n".format(unassigned))
    handle.write("__bim_rows__\\t{}\\n".format(len(variants)))

# Written here rather than captured by an `eval` output, which Nextflow allows only on a Bash script.
with open("versions.yml", "w", newline="\\n") as handle:
    handle.write('"{}":\\n    python: {}\\n'.format(process_name, sys.version.split()[0]))
