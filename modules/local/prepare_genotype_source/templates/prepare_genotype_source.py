#!/usr/bin/env python3
"""Digest a supplied genotype bundle where its bytes live.

The identity a multi-gigabyte bundle gets must not be computed in the Nextflow JVM on the head node, and it
must survive `-resume`, so it is computed here in a task. Two artifacts come out: the per-member table, which
is provenance, and the single bundle digest, which is the identity everything downstream keys on.
"""

import hashlib
import json
import os
import sys

MEMBERS = json.loads(${members_literal})
PREFIX = ${prefix_literal}
PROCESS_NAME = ${task_process_literal}

# A member's role is its extension with any compression suffix removed: `pgen`, `psam`, `pvar`, `bed`, `bim`,
# or `fam`. Deriving it here rather than passing it in keeps a copy of the manifest's genotype-group
# contract out of this module; the caller owns the order the members arrive in.
COMPRESSION_SUFFIXES = (".gz", ".bgz", ".zst")
CHUNK_BYTES = 1 << 20


def role_of(name):
    stem = name
    for suffix in COMPRESSION_SUFFIXES:
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return stem.rsplit(".", 1)[-1]


def digest_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


rows = []
for member in MEMBERS:
    name = os.path.basename(member)
    rows.append((role_of(name), name, str(os.path.getsize(member)), digest_of(member)))

with open("{}.genotype_source.tsv".format(PREFIX), "w", encoding="utf-8") as handle:
    for row in rows:
        handle.write("\\t".join(row) + "\\n")

# The bundle identity folds only the roles and their digests, in the order the caller supplied. File names and
# byte sizes are deliberately excluded: renaming or moving a byte-identical bundle must not change the
# identity that published artifact keys are built from.
identity = hashlib.sha256()
for role, _name, _size, sha256 in rows:
    identity.update("{}\\t{}\\n".format(role, sha256).encode("utf-8"))

with open("{}.genotype_source.sha256".format(PREFIX), "w", encoding="utf-8") as handle:
    handle.write(identity.hexdigest() + "\\n")

with open("versions.yml", "w", newline="") as handle:
    handle.write('"{}":\\n'.format(PROCESS_NAME))
    handle.write("    python: {}\\n".format(sys.version.split()[0]))
