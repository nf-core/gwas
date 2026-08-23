#!/usr/bin/env python3

import csv
import json
import math
import re
import sys
from pathlib import Path


META = json.loads(${meta_literal})
HERS = Path(${hers_literal})
EXTRA = Path(${extra_literal})
OVERLAP = Path(${overlap_literal})
LOG = Path(${log_literal})
PREPARATION = Path(${preparation_literal})
HERS_LIABILITY = Path(${hers_liability_literal}) if ${hers_liability_literal} else None
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
            values[fields[0]] = fields[1]
    return values


def her_all(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = {row["Component"]: row for row in csv.DictReader(handle, delimiter=" ", skipinitialspace=True)}
    if "Her_All" not in rows:
        raise ValueError(f"mandatory Her_All row is absent from {path}")
    return number(rows["Her_All"].get("Heritability")), number(rows["Her_All"].get("SE"))


observed_estimate, observed_se = her_all(HERS)
extra = key_values(EXTRA)
overlap = key_values(OVERLAP)
log_text = LOG.read_text(encoding="utf-8")
if extra.get("Converged") != "YES":
    raise ValueError("LDAK SumHer did not report Converged YES")

warnings = [line.strip() for line in log_text.splitlines() if line.strip().startswith("Warning,")]
classification = "estimable"
if observed_estimate is None or observed_se is None:
    classification = "completed_nonestimable"
elif not 0 <= observed_estimate <= 1 or warnings:
    classification = "estimable_with_warning"

rows = [["observed", observed_estimate, observed_se]]
if HERS_LIABILITY and HERS_LIABILITY.exists():
    liability_estimate, liability_se = her_all(HERS_LIABILITY)
    rows.append(["liability", liability_estimate, liability_se])
    if liability_estimate is None or liability_se is None:
        classification = "completed_nonestimable"
    elif classification == "estimable" and not 0 <= liability_estimate <= 1:
        classification = "estimable_with_warning"

def rendered(value):
    return "NA" if value is None else format(value, ".17g")


with Path("heritability.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t", lineterminator="\\n")
    writer.writerow(["summary_statistics_id", "request_id", "method", "trait_id", "trait_type", "scale", "estimate", "standard_error", "classification"])
    for scale, estimate, standard_error in rows:
        writer.writerow([META["summary_statistics_id"], META["request_id"], META["method"], META["trait_id"], META["trait_type"], scale, rendered(estimate), rendered(standard_error), classification])

effective_args = META.get("effective_native_args", [])
option_names = [token.split("=", 1)[0] for token in effective_args if isinstance(token, str) and token.startswith("--")]
large_effect_match = re.search(r"There are (\\d+) predictors that explain at least .*?; these will be excluded", log_text)
overlap_proportion = number(overlap.get("Overlap_Proportion"))
preparation = json.loads(PREPARATION.read_text(encoding="utf-8"))
if "--check-sums" in option_names and "NO" in effective_args and overlap_proportion is not None and overlap_proportion < 0.8:
    warnings.append("HIGH: summary-statistics coverage is below LDAK's approximate 80% guidance under an explicit incomplete-summary override")

diagnostics = {
    "classification": classification,
    "converged": extra.get("Converged"),
    "tagging_file_predictors": overlap.get("Tagging_File_Predictors"),
    "summary_statistic_predictors": overlap.get("Summary_Statistic_Predictors"),
    "overlap_proportion": overlap.get("Overlap_Proportion"),
    "predictors_used": extra.get("Num_Datapoints"),
    "large_effect_predictors_excluded": int(large_effect_match.group(1)) if large_effect_match else None,
    "strand_ambiguous_rows_present_before_native_filtering": preparation.get("strand_ambiguous_rows_present"),
    "weighted_sample_size": extra.get("Weighted_Sample_Size"),
    "effective_size": extra.get("Effective_Size"),
    "weighted_mean_chisq": extra.get("Weighted_Mean_Chisq"),
    "weighted_gif": extra.get("Weighted_GIF"),
    "null_log_likelihood": extra.get("Null_logl"),
    "alternative_log_likelihood": extra.get("Alt_logl"),
}
with Path("diagnostics.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t", lineterminator="\\n")
    writer.writerow(["request_id", "metric", "value"])
    for metric, value in diagnostics.items():
        writer.writerow([META["request_id"], metric, "NA" if value is None else value])

provenance = {
    "schema_version": "1.0",
    "request_id": META["request_id"],
    "primary_request_id": META.get("primary_request_id"),
    "request_name": META.get("request_name"),
    "method": META["method"],
    "summary_statistics_id": META["summary_statistics_id"],
    "trait": {"trait_id": META["trait_id"], "trait_type": META["trait_type"]},
    "declared_prevalence": {"population": META.get("population_prevalence"), "sample": META.get("sample_prevalence")},
    "reference": META.get("reference_metadata"),
    "reference_validation": META.get("reference_validation"),
    "native_args": META.get("native_args", []),
    "effective_native_args": effective_args,
    "scientific_policy": {
        "large_effect": "truncate" if "--truncate" in option_names else "cutoff",
        "ambiguous_variants": "retained_by_explicit_override" if "--allow-ambiguous" in option_names else "excluded_by_native_default",
        "complete_summary_check": "disabled_by_explicit_override" if "--check-sums" in option_names and "NO" in effective_args else "enabled_by_native_default",
        "scope": "total_common_snp_heritability",
    },
    "preparation": preparation,
    "diagnostics": diagnostics,
    "warnings": warnings,
    "classification": classification,
    "native_artifacts": [
        "native.hers",
        "native.cats",
        "native.share",
        "native.enrich",
        "native.extra",
        "native.cross",
        "native.taus",
        "native.labels",
        "native.progress",
        "native.overlap",
        "native.log",
    ] + (["native.hers.liab", "native.cats.liab", "native.factor"] if HERS_LIABILITY and HERS_LIABILITY.exists() else []),
    "execution": {
        "native": {"owner": "LDAK_SUMHER", "declared_container": META.get("native_runtime"), "tool": "LDAK", "version": "6.3"},
        "normalization": {"owner": PROCESS, "container": CONTAINER, "tool": "Python", "version": sys.version.split()[0]},
    },
}
Path("provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\\n", encoding="utf-8")
Path("versions.yml").write_text(f'"{PROCESS}":\\n    python: {sys.version.split()[0]}\\n', encoding="utf-8")
