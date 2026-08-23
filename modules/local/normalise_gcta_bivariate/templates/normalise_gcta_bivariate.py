#!/usr/bin/env python3
"""Normalize native dense GCTA bivariate REML outputs without choosing a preferred result."""

import json
import math
import re
import sys

META = json.loads(${meta_literal})
HSQ = ${hsq_literal}
GCTA_LOG = ${gcta_log_literal}
PAIR_LOG = ${pair_log_literal}
PROCESS_NAME = ${task_process_literal}
CONTAINER = ${task_container_literal}
METHOD = "gcta_bivariate_reml"


def fail(message):
    sys.exit("[nf-core/gwas] ERROR: pair request '{}': {}".format(META["request_id"], message))


def parse_number(text):
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def display(value):
    return "NA" if value is None else format(value, ".15g")


def parse_hsq(path):
    components = {}
    with open(path) as handle:
        lines = [line.strip() for line in handle if line.strip()]
    for line in lines:
        fields = line.split()
        if fields[:3] == ["Source", "Variance", "SE"]:
            continue
        if len(fields) not in (2, 3):
            continue
        source = fields[0]
        if source in components:
            fail("native result '{}' repeats component '{}'".format(path, source))
        estimate_text = fields[1]
        se_text = fields[2] if len(fields) == 3 else None
        components[source] = {
            "estimate": parse_number(estimate_text),
            "standard_error": parse_number(se_text) if se_text is not None else None,
            "native_estimate": estimate_text,
            "native_standard_error": se_text,
        }
    mandatory = ["V(G)_tr1", "V(G)_tr2", "C(G)_tr12", "V(e)_tr1", "V(e)_tr2", "Vp_tr1", "Vp_tr2", "V(G)/Vp_tr1", "V(G)/Vp_tr2", "rG", "logL", "n"]
    missing = [source for source in mandatory if source not in components]
    if missing:
        fail("native result '{}' is missing mandatory components {}".format(path, ", ".join(missing)))
    return components


def parse_pair_log(path):
    values = {}
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\\r\\n").split("\\t")
            if len(fields) != 2 or not fields[0]:
                fail("pair diagnostic '{}' line {} must contain one metric and one value".format(path, line_number))
            values[fields[0]] = fields[1]
    required = ["left_samples", "right_samples", "endpoint_overlap_samples", "union_samples", "left_nonmissing", "right_nonmissing", "both_nonmissing", "quantitative_covariate_samples", "categorical_covariate_samples"]
    missing = [name for name in required if name not in values]
    if missing:
        fail("pair diagnostic '{}' is missing metrics {}".format(path, ", ".join(missing)))
    return values


with open(GCTA_LOG) as handle:
    log_text = handle.read()
components = parse_hsq(HSQ)
pair_diagnostics = parse_pair_log(PAIR_LOG)

converged = bool(re.search(r"Log-likelihood ratio converged", log_text, flags=re.IGNORECASE))
common_match = re.search(r"(\\d+) individuals are in common", log_text)
nonmissing_match = re.search(r"(\\d+) non-missing phenotypes for trait #1 and (\\d+) for trait #2", log_text)
version_match = re.search(r"version v([0-9.]+)", log_text)

if not converged:
    fail("native log '{}' does not confirm REML convergence".format(GCTA_LOG))

warning_patterns = [
    ("bended_information_matrix", r"matrix is bended"),
    ("constrained_component", r"component\\(s\\) constrained|constrain the .* component"),
    ("unreliable_standard_error", r"SE is unreliable"),
]
warnings = [name for name, pattern in warning_patterns if re.search(pattern, log_text, flags=re.IGNORECASE)]
rg = components["rG"]["estimate"]
rg_se = components["rG"]["standard_error"]
left_h2 = components["V(G)/Vp_tr1"]["estimate"]
right_h2 = components["V(G)/Vp_tr2"]["estimate"]
if rg is not None and abs(rg) > 1:
    warnings.append("genetic_correlation_outside_unit_interval")
elif rg is not None and abs(rg) == 1:
    warnings.append("genetic_correlation_at_boundary")

for side, h2 in [("left", left_h2), ("right", right_h2)]:
    if h2 is not None and (h2 < 0 or h2 > 1):
        warnings.append(side + "_heritability_outside_unit_interval")
    elif h2 is not None and h2 in (0, 1):
        warnings.append(side + "_heritability_at_boundary")

for source in ["V(G)_tr1", "V(G)_tr2", "V(e)_tr1", "V(e)_tr2", "Vp_tr1", "Vp_tr2"]:
    if components[source]["estimate"] is not None and components[source]["estimate"] < 0:
        warnings.append("negative_variance_component:" + source)

primary_components = [components[name] for name in ["V(G)/Vp_tr1", "V(G)/Vp_tr2", "C(G)_tr12", "rG"]]
estimable = all(component["estimate"] is not None and component["standard_error"] is not None for component in primary_components)
if not converged or not estimable:
    classification = "completed_nonestimable"
elif warnings:
    classification = "estimable_with_warning"
else:
    classification = "estimable"


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        handle.write("\\t".join(header) + "\\n")
        for row in rows:
            handle.write("\\t".join(str(value) for value in row) + "\\n")


common = [META["relationship_id"], META["request_id"], METHOD]
heritability_rows = []
for side, suffix in [("left", "tr1"), ("right", "tr2")]:
    endpoint = {
        "analysis_id": META[side + "_analysis_id"],
        "trait_id": META[side + "_trait_id"],
        "trait_type": META[side + "_trait_type"],
    }
    observed = components["V(G)/Vp_" + suffix]
    heritability_rows.append(common + [side, endpoint["analysis_id"], endpoint["trait_id"], endpoint["trait_type"], "observed", display(observed["estimate"]), display(observed["standard_error"]), classification])
    liability_key = "V(G)/Vp_{}_L".format(suffix)
    if liability_key in components:
        liability = components[liability_key]
        heritability_rows.append(common + [side, endpoint["analysis_id"], endpoint["trait_id"], endpoint["trait_type"], "liability", display(liability["estimate"]), display(liability["standard_error"]), classification])
write_tsv(
    "heritability.tsv",
    ["relationship_id", "request_id", "method", "endpoint", "analysis_id", "trait_id", "trait_type", "scale", "estimate", "standard_error", "classification"],
    heritability_rows,
)

write_tsv(
    "genetic_correlation.tsv",
    ["relationship_id", "request_id", "method", "left_analysis_id", "right_analysis_id", "left_trait_id", "right_trait_id", "estimate", "standard_error", "classification"],
    [common + [META["left_analysis_id"], META["right_analysis_id"], META["left_trait_id"], META["right_trait_id"], display(rg), display(rg_se), classification]],
)
covariance = components["C(G)_tr12"]
write_tsv(
    "genetic_covariance.tsv",
    ["relationship_id", "request_id", "method", "left_analysis_id", "right_analysis_id", "left_trait_id", "right_trait_id", "scale", "estimate", "standard_error", "classification"],
    [common + [META["left_analysis_id"], META["right_analysis_id"], META["left_trait_id"], META["right_trait_id"], "observed", display(covariance["estimate"]), display(covariance["standard_error"]), classification]],
)

diagnostics = dict(pair_diagnostics)
residual_covariance_present = "C(e)_tr12" in components
if residual_covariance_present:
    residual_covariance_status = "retained"
elif "--reml-bivar-nocove" in META.get("native_args", []):
    residual_covariance_status = "dropped_by_native_option"
else:
    residual_covariance_status = "dropped_by_native_overlap_rule"
diagnostics.update({
    "classification": classification,
    "native_converged": str(converged).lower(),
    "native_common_samples": common_match.group(1) if common_match else "NA",
    "native_trait1_nonmissing": nonmissing_match.group(1) if nonmissing_match else "NA",
    "native_trait2_nonmissing": nonmissing_match.group(2) if nonmissing_match else "NA",
    "native_observations": components["n"]["native_estimate"],
    "native_log_likelihood": components["logL"]["native_estimate"],
    "residual_covariance": components.get("C(e)_tr12", {}).get("native_estimate", "NA"),
    "residual_covariance_status": residual_covariance_status,
    "warnings": ",".join(sorted(set(warnings))) if warnings else "none",
})
write_tsv(
    "diagnostics.tsv",
    ["relationship_id", "request_id", "metric", "value"],
    [[META["relationship_id"], META["request_id"], name, value] for name, value in sorted(diagnostics.items())],
)

provenance = {
    "schema_version": "1.0",
    "relationship_id": META["relationship_id"],
    "request_id": META["request_id"],
    "method": METHOD,
    "orientation": {
        "left": {"analysis_id": META["left_analysis_id"], "trait_id": META["left_trait_id"], "trait_type": META["left_trait_type"]},
        "right": {"analysis_id": META["right_analysis_id"], "trait_id": META["right_trait_id"], "trait_type": META["right_trait_type"]},
    },
    "cohort_id": META["cohort"],
    "matrix": {
        "kind": META["matrix_kind"],
        "key": META["matrix_key"],
        "basename": META["matrix_basename"],
        "settings": META["matrix_settings"],
    },
    "effective_prevalence": META["reml_bivar_prevalence"],
    "declared_prevalence": {
        "left": META["left_population_prevalence"],
        "right": META["right_population_prevalence"],
    },
    "native_args": META["native_args"],
    "classification": classification,
    "warnings": sorted(set(warnings)),
    "diagnostics": diagnostics,
    "native_components": components,
    "execution": {"process": PROCESS_NAME, "container": CONTAINER},
    "tool": {"name": "gcta", "version": version_match.group(1) if version_match else None},
    "artifacts": {
        "native": ["native.hsq", "native.log"],
        "normalized": ["heritability.tsv", "genetic_correlation.tsv", "genetic_covariance.tsv", "diagnostics.tsv", "provenance.json"],
    },
}
with open("provenance.json", "w", newline="") as handle:
    json.dump(provenance, handle, indent=2, sort_keys=True, allow_nan=False)
    handle.write("\\n")

with open("versions.yml", "w", newline="") as handle:
    handle.write('"{}":\\n'.format(PROCESS_NAME))
    handle.write("    python: {}\\n".format(sys.version.split()[0]))
