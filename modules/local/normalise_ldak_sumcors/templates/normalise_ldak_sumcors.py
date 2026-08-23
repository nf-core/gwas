#!/usr/bin/env python3

import csv
import json
import math
import re
import sys
from pathlib import Path


META = json.loads(${meta_literal})
CORRELATIONS = Path(${correlations_literal})
CORRELATIONS_FULL = Path(${correlations_full_literal})
OVERLAP = Path(${overlap_literal})
LOG = Path(${log_literal})
LEFT_PREPARATION = Path(${left_preparation_literal})
RIGHT_PREPARATION = Path(${right_preparation_literal})
CORRELATIONS_LIABILITY = Path(${correlations_liability_literal}) if ${correlations_liability_literal} else None
PROCESS = ${task_process_literal}
CONTAINER = ${task_container_literal}


def number(value):
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def key_values(path):
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 2:
            values[fields[0]] = fields[1:]
    return values


native = key_values(CORRELATIONS)
mandatory = ["Her1_All", "Her2_All", "Coher_All", "Cor_All", "Intercept1", "Intercept2", "Overlap"]
missing = [component for component in mandatory if component not in native]
if missing:
    raise ValueError(f"mandatory SumCors components are absent: {', '.join(missing)}")

observed = {component: [number(value) for value in native[component][:2]] for component in mandatory}
log_text = LOG.read_text(encoding="utf-8")
warnings = [line.strip() for line in log_text.splitlines() if line.strip().startswith("Warning,")]
primary = [observed[component][index] for component in ["Her1_All", "Her2_All", "Coher_All", "Cor_All"] for index in [0, 1]]
classification = "completed_nonestimable" if any(value is None for value in primary) else "estimable"
if classification == "estimable" and (not 0 <= observed["Her1_All"][0] <= 1 or not 0 <= observed["Her2_All"][0] <= 1 or not -1 <= observed["Cor_All"][0] <= 1 or warnings):
    classification = "estimable_with_warning"

liability = None
if CORRELATIONS_LIABILITY and CORRELATIONS_LIABILITY.exists():
    liability_native = key_values(CORRELATIONS_LIABILITY)
    liability = {component: [number(value) for value in liability_native.get(component, [])[:2]] for component in ["Her1_All", "Her2_All", "Coher_All", "Cor_All"]}
    if any(len(values) < 2 or value is None for values in liability.values() for value in values):
        classification = "completed_nonestimable"


def rendered(value):
    return "NA" if value is None else format(value, ".17g")


with Path("heritability.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t", lineterminator="\\n")
    writer.writerow(["relationship_id", "request_id", "method", "endpoint", "summary_statistics_id", "trait_id", "trait_type", "scale", "estimate", "standard_error", "classification"])
    for scale, values in [("observed", observed)] + ([('liability', liability)] if liability else []):
        for endpoint, component in [("left", "Her1_All"), ("right", "Her2_All")]:
            prefix = "left" if endpoint == "left" else "right"
            writer.writerow([META["relationship_id"], META["request_id"], META["method"], endpoint, META[f"{prefix}_summary_statistics_id"], META[f"{prefix}_trait_id"], META[f"{prefix}_trait_type"], scale, rendered(values[component][0]), rendered(values[component][1]), classification])

with Path("genetic_covariance.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t", lineterminator="\\n")
    writer.writerow(["relationship_id", "request_id", "method", "left_summary_statistics_id", "right_summary_statistics_id", "left_trait_id", "right_trait_id", "scale", "estimate", "standard_error", "classification"])
    for scale, values in [("observed", observed)] + ([('liability', liability)] if liability else []):
        writer.writerow([META["relationship_id"], META["request_id"], META["method"], META["left_summary_statistics_id"], META["right_summary_statistics_id"], META["left_trait_id"], META["right_trait_id"], scale, rendered(values["Coher_All"][0]), rendered(values["Coher_All"][1]), classification])

with Path("genetic_correlation.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t", lineterminator="\\n")
    writer.writerow(["relationship_id", "request_id", "method", "left_summary_statistics_id", "right_summary_statistics_id", "left_trait_id", "right_trait_id", "estimate", "standard_error", "classification"])
    writer.writerow([META["relationship_id"], META["request_id"], META["method"], META["left_summary_statistics_id"], META["right_summary_statistics_id"], META["left_trait_id"], META["right_trait_id"], rendered(observed["Cor_All"][0]), rendered(observed["Cor_All"][1]), classification])

overlap = key_values(OVERLAP)
left_preparation = json.loads(LEFT_PREPARATION.read_text(encoding="utf-8"))
right_preparation = json.loads(RIGHT_PREPARATION.read_text(encoding="utf-8"))
effective_args = META.get("effective_native_args", [])
option_names = [token.split("=", 1)[0] for token in effective_args if isinstance(token, str) and token.startswith("--")]
left_excluded_match = re.search(r"There are (\\d+) predictors that explain at least .*? for Trait 1; these will be excluded", log_text)
right_excluded_match = re.search(r"There are (\\d+) predictors that explain at least .*? for Trait 2; these will be excluded", log_text)
predictors_used_match = re.search(r"usually tens of thousands of predictors \\(not (\\d+)\\)", log_text)
overlap_proportion = number(overlap.get("Overlap_Proportion", [None])[0])
if "--check-sums" in option_names and "NO" in effective_args and overlap_proportion is not None and overlap_proportion < 0.8:
    warnings.append("HIGH: summary-statistics coverage is below LDAK's approximate 80% guidance under an explicit incomplete-summary override")
diagnostics = {
    "classification": classification,
    "tagging_file_predictors": overlap.get("Tagging_File_Predictors", [None])[0],
    "summary_statistic_predictors": overlap.get("Summary_Statistic_Predictors", [None])[0],
    "overlap_proportion": overlap.get("Overlap_Proportion", [None])[0],
    "predictors_used": int(predictors_used_match.group(1)) if predictors_used_match else None,
    "left_large_effect_predictors_excluded": int(left_excluded_match.group(1)) if left_excluded_match else None,
    "right_large_effect_predictors_excluded": int(right_excluded_match.group(1)) if right_excluded_match else None,
    "left_strand_ambiguous_rows_present_before_native_filtering": left_preparation.get("strand_ambiguous_rows_present"),
    "right_strand_ambiguous_rows_present_before_native_filtering": right_preparation.get("strand_ambiguous_rows_present"),
    "left_intercept": rendered(observed["Intercept1"][0]),
    "left_intercept_se": rendered(observed["Intercept1"][1]),
    "right_intercept": rendered(observed["Intercept2"][0]),
    "right_intercept_se": rendered(observed["Intercept2"][1]),
    "native_overlap_estimate": rendered(observed["Overlap"][0]),
    "native_overlap_standard_error": rendered(observed["Overlap"][1]),
}
with Path("diagnostics.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t", lineterminator="\\n")
    writer.writerow(["request_id", "metric", "value"])
    for metric, value in diagnostics.items():
        writer.writerow([META["request_id"], metric, "NA" if value is None else value])

provenance = {
    "schema_version": "1.0",
    "relationship_id": META["relationship_id"],
    "request_id": META["request_id"],
    "primary_request_id": META.get("primary_request_id"),
    "request_name": META.get("request_name"),
    "method": META["method"],
    "orientation": {
        "left": {"summary_statistics_id": META["left_summary_statistics_id"], "trait_id": META["left_trait_id"], "trait_type": META["left_trait_type"]},
        "right": {"summary_statistics_id": META["right_summary_statistics_id"], "trait_id": META["right_trait_id"], "trait_type": META["right_trait_type"]},
    },
    "declared_prevalence": {
        "left": {"population": META.get("left_population_prevalence"), "sample": META.get("left_sample_prevalence")},
        "right": {"population": META.get("right_population_prevalence"), "sample": META.get("right_sample_prevalence")},
    },
    "reference": META.get("reference_metadata"),
    "reference_validation": META.get("reference_validation"),
    "native_args": META.get("native_args", []),
    "effective_native_args": effective_args,
    "scientific_policy": {
        "large_effect": "truncate" if "--truncate" in option_names else "cutoff",
        "ambiguous_variants": "retained_by_explicit_override" if "--allow-ambiguous" in option_names else "excluded_by_native_default",
        "complete_summary_check": "disabled_by_explicit_override" if "--check-sums" in option_names and "NO" in effective_args else "enabled_by_native_default",
        "sample_overlap": "constrained_by_explicit_override" if "--overlapping-samples" in option_names else "estimated_by_native_default",
        "scope": "total_common_snp_genetic_correlation",
    },
    "preparation": {
        "left": left_preparation,
        "right": right_preparation,
    },
    "diagnostics": diagnostics,
    "warnings": warnings,
    "classification": classification,
    "native_artifacts": [
        "native.cors",
        "native.cors.full",
        "native.labels",
        "native.progress",
        "native.overlap",
        "native.log",
    ] + (["native.cors.liab"] if CORRELATIONS_LIABILITY and CORRELATIONS_LIABILITY.exists() else []),
    "execution": {
        "native": {"owner": "LDAK_SUMCORS", "declared_container": META.get("native_runtime"), "tool": "LDAK", "version": "6.3"},
        "normalization": {"owner": PROCESS, "container": CONTAINER, "tool": "Python", "version": sys.version.split()[0]},
    },
}
Path("provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\\n", encoding="utf-8")
Path("versions.yml").write_text(f'"{PROCESS}":\\n    python: {sys.version.split()[0]}\\n', encoding="utf-8")
