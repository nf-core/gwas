#!/usr/bin/env python3
"""Normalize native CBIIT LDSC logs without altering or ranking estimates."""

import json
import math
import re
import sys

META = json.loads(${meta_literal})
OBSERVED_LOG = ${observed_log_literal}
LIABILITY_LOG = ${liability_log_literal}
MUNGING_LOGS = ${munging_logs_literal}
PROCESS_NAME = ${task_process_literal}
CONTAINER = ${task_container_literal}
LDSC_CONTAINER = ${ldsc_container_literal}
METHOD = META["method"]


def fail(message):
    sys.exit("[nf-core/gwas] ERROR: request '{}': {}".format(META["request_id"], message))


def parse_number(text):
    if text is None:
        return None
    cleaned = text.strip().rstrip(".")
    if cleaned.lower() in {"na", "nan", "inf", "+inf", "-inf"}:
        return None
    try:
        value = float(cleaned)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def display(value):
    return "NA" if value is None else format(value, ".15g")


def payload_for(text, label):
    prefix = label.lower() + ":"
    for line in text.splitlines():
        if line.lower().startswith(prefix):
            return line.split(":", 1)[1].strip()
    return None


def find_value(text, label, required=True):
    payload = payload_for(text, label)
    if payload is None:
        if required:
            fail("native log is missing '{}'".format(label))
        return {"estimate": None, "native": None}
    native = payload.split()[0] if payload else None
    return {"estimate": parse_number(native), "native": native}


def find_value_se(text, label, required=True):
    payload = payload_for(text, label)
    left = payload.rfind("(") if payload is not None else -1
    right = payload.rfind(")") if payload is not None else -1
    if payload is None or left < 1 or right <= left:
        if required:
            fail("native log is missing '{}' with its standard error".format(label))
        return {"estimate": None, "standard_error": None, "native_estimate": None, "native_standard_error": None}
    native_estimate = payload[:left].strip().split()[0]
    native_standard_error = payload[left + 1 : right].strip()
    return {
        "estimate": parse_number(native_estimate),
        "standard_error": parse_number(native_standard_error),
        "native_estimate": native_estimate,
        "native_standard_error": native_standard_error,
    }


def log_scale(text, estimand):
    for scale in ("observed", "liability"):
        label = "Total {} scale {}".format(scale.title(), estimand)
        payload = payload_for(text, label)
        if payload is not None:
            parsed = find_value_se(text, label)
            return scale, parsed
    fail("native log is missing total {} and its standard error".format(estimand))


def warnings_from(text):
    return sorted(set(line.split(":", 1)[1].strip() for line in text.splitlines() if line.lower().startswith("warning:")))


def version_from(text):
    for line in text.splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[0] == "*" and fields[1] == "Version":
            return fields[2]
    return None


def regression_snps(text, pairwise=False):
    if pairwise:
        valid = re.findall(r"^([0-9]+) SNPs with valid alleles[.]", text, flags=re.MULTILINE)
        if valid:
            return int(valid[-1])
    merged = re.findall(r"^After merging with regression SNP LD, ([0-9]+) SNPs remain[.]", text, flags=re.MULTILINE)
    if not merged:
        fail("native log does not report the regression-SNP count")
    return int(merged[-1])


def parse_munging(path):
    with open(path) as handle:
        text = handle.read()
    read = re.search(r"Read summary statistics for ([0-9]+) SNPs", text)
    retained = re.findall(r"([0-9]+) SNPs remain", text)
    return {
        "file": path,
        "input_snps": int(read.group(1)) if read else None,
        "retained_snps": int(retained[-1]) if retained else None,
        "warnings": warnings_from(text),
    }


def classify(values, warnings, bounded):
    finite = all(value["estimate"] is not None and value.get("standard_error") is not None for value in values)
    if not finite:
        return "completed_nonestimable", warnings
    warning_codes = list(warnings)
    for label, value, lower, upper in bounded:
        if value is not None and (value < lower or value > upper):
            warning_codes.append("{}_outside_{}_{}".format(label, lower, upper))
        elif value is not None and value in (lower, upper):
            warning_codes.append("{}_at_boundary".format(label))
    return ("estimable_with_warning" if warning_codes else "estimable"), sorted(set(warning_codes))


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        handle.write(chr(9).join(header) + chr(10))
        for row in rows:
            handle.write(chr(9).join(str(value) for value in row) + chr(10))


def diagnostic_rows(classification, warnings, scales, native):
    diagnostics = {
        "classification": classification,
        "native_scales": ",".join(scales),
        "regression_snps": native["regression_snps"],
        "warnings": " | ".join(warnings) if warnings else "none",
    }
    if METHOD == "ldsc_h2":
        diagnostics.update(
            {
                "lambda_gc": display(native["lambda_gc"]["estimate"]),
                "mean_chi_squared": display(native["mean_chi_squared"]["estimate"]),
                "intercept": display(native["intercept"]["estimate"]),
                "intercept_standard_error": display(native["intercept"]["standard_error"]),
                "ratio": display(native["ratio"]["estimate"]),
                "ratio_standard_error": display(native["ratio"]["standard_error"]),
            }
        )
    else:
        diagnostics.update(
            {
                "left_intercept": display(native["left_intercept"]["estimate"]),
                "left_intercept_standard_error": display(native["left_intercept"]["standard_error"]),
                "right_intercept": display(native["right_intercept"]["estimate"]),
                "right_intercept_standard_error": display(native["right_intercept"]["standard_error"]),
                "cross_trait_intercept": display(native["cross_trait_intercept"]["estimate"]),
                "cross_trait_intercept_standard_error": display(native["cross_trait_intercept"]["standard_error"]),
                "mean_z1_z2": display(native["mean_z1_z2"]["estimate"]),
            }
        )
    relationship = META.get("relationship_id") or "NA"
    return [[relationship, META["request_id"], key, value] for key, value in sorted(diagnostics.items())]


with open(OBSERVED_LOG) as handle:
    observed_text = handle.read()
liability_text = None
if LIABILITY_LOG:
    with open(LIABILITY_LOG) as handle:
        liability_text = handle.read()

texts = [(OBSERVED_LOG, observed_text)] + ([(LIABILITY_LOG, liability_text)] if liability_text is not None else [])
native_version = version_from(observed_text)
if not native_version:
    fail("observed-scale native log does not report the LDSC version")
if liability_text is not None and version_from(liability_text) != native_version:
    fail("observed- and liability-scale logs report different LDSC versions")

munging = [parse_munging(path) for path in MUNGING_LOGS]
all_native_warnings = sorted(set(warning for _path, text in texts for warning in warnings_from(text)))

if METHOD == "ldsc_h2":
    estimates = []
    native_by_scale = {}
    for path, text in texts:
        scale, h2 = log_scale(text, "h2")
        if scale in native_by_scale:
            fail("native inputs repeat '{}' scale H2".format(scale))
        native_by_scale[scale] = {"h2": h2, "artifact": "native.{}.log".format(scale)}
        estimates.append(h2)
    if "observed" not in native_by_scale:
        fail("the primary LDSC H2 invocation did not emit observed-scale H2")
    primary_text = observed_text
    native = {
        "regression_snps": regression_snps(primary_text),
        "lambda_gc": find_value(primary_text, "Lambda GC"),
        "mean_chi_squared": find_value(primary_text, "Mean Chi^2"),
        "intercept": find_value_se(primary_text, "Intercept"),
        "ratio": find_value_se(primary_text, "Ratio", required=False),
        "estimates": native_by_scale,
    }
    classification, warning_codes = classify(
        estimates,
        all_native_warnings,
        [("heritability_{}".format(scale), values["h2"]["estimate"], 0, 1) for scale, values in native_by_scale.items()],
    )
    rows = []
    for scale, values in native_by_scale.items():
        h2 = values["h2"]
        rows.append(
            [
                "NA",
                META["request_id"],
                METHOD,
                "unary",
                META["summary_statistics_id"],
                META["trait_id"],
                META["trait_type"],
                scale,
                display(h2["estimate"]),
                display(h2["standard_error"]),
                classification,
                "requests/{}/{}/{}".format(METHOD, META["request_id"], values["artifact"]),
            ]
        )
    write_tsv(
        "heritability.tsv",
        ["relationship_id", "request_id", "method", "endpoint", "summary_statistics_id", "trait_id", "trait_type", "scale", "estimate", "standard_error", "classification", "native_artifact"],
        rows,
    )
else:
    if METHOD != "ldsc_rg":
        fail("unsupported LDSC normalization method '{}'".format(METHOD))

    def section(text, start, end):
        lower_text = text.lower()
        start_at = lower_text.find(start.lower())
        end_at = lower_text.find(end.lower(), start_at + len(start)) if start_at >= 0 else -1
        if start_at < 0 or end_at < 0:
            fail("native RG log is missing the '{}' section".format(start))
        return text[start_at + len(start) : end_at]

    native_by_scale = {}
    primary_correlation = None
    primary_diagnostics = None
    estimands = []
    for path, text in texts:
        left_text = section(text, r"Heritability of phenotype 1", r"Heritability of phenotype 2/2")
        right_text = section(text, r"Heritability of phenotype 2/2", r"Genetic Covariance")
        covariance_text = section(text, r"Genetic Covariance", r"Genetic Correlation")
        correlation_text = section(text, r"Genetic Correlation", r"Summary of Genetic Correlation Results")
        left_scale, left_h2 = log_scale(left_text, "h2")
        right_scale, right_h2 = log_scale(right_text, "h2")
        covariance_scale, covariance = log_scale(covariance_text, "gencov")
        if len({left_scale, right_scale, covariance_scale}) != 1:
            fail("native RG log mixes estimand scales within one invocation")
        scale = left_scale
        if scale in native_by_scale:
            fail("native inputs repeat '{}' scale RG".format(scale))
        correlation = find_value_se(correlation_text, "Genetic Correlation")
        z_score = find_value(correlation_text, "Z-score")
        p_value = find_value(correlation_text, "P")
        native_by_scale[scale] = {
            "left_h2": left_h2,
            "right_h2": right_h2,
            "covariance": covariance,
            "correlation": correlation,
            "z_score": z_score,
            "p_value": p_value,
            "artifact": "native.{}.log".format(scale),
        }
        estimands.extend([left_h2, right_h2, covariance, correlation])
        if scale == "observed":
            primary_correlation = native_by_scale[scale]
            primary_diagnostics = {
                "regression_snps": regression_snps(text, pairwise=True),
                "left_intercept": find_value_se(left_text, "Intercept"),
                "right_intercept": find_value_se(right_text, "Intercept"),
                "cross_trait_intercept": find_value_se(covariance_text, "Intercept"),
                "mean_z1_z2": find_value(covariance_text, "Mean z1*z2"),
            }
    if primary_correlation is None:
        fail("the primary LDSC RG invocation did not emit observed-scale estimands")
    classification, warning_codes = classify(
        estimands,
        all_native_warnings,
        [("left_heritability_{}".format(scale), values["left_h2"]["estimate"], 0, 1) for scale, values in native_by_scale.items()]
        + [("right_heritability_{}".format(scale), values["right_h2"]["estimate"], 0, 1) for scale, values in native_by_scale.items()]
        + [("genetic_correlation", primary_correlation["correlation"]["estimate"], -1, 1)],
    )
    heritability_rows = []
    covariance_rows = []
    for scale, values in native_by_scale.items():
        artifact = "requests/{}/{}/{}".format(METHOD, META["request_id"], values["artifact"])
        for endpoint, summary_key, trait_key, type_key, value_key in [
            ("left", "left_summary_statistics_id", "left_trait_id", "left_trait_type", "left_h2"),
            ("right", "right_summary_statistics_id", "right_trait_id", "right_trait_type", "right_h2"),
        ]:
            value = values[value_key]
            heritability_rows.append(
                [META["relationship_id"], META["request_id"], METHOD, endpoint, META[summary_key], META[trait_key], META[type_key], scale, display(value["estimate"]), display(value["standard_error"]), classification, artifact]
            )
        covariance = values["covariance"]
        covariance_rows.append(
            [META["relationship_id"], META["request_id"], METHOD, META["left_summary_statistics_id"], META["right_summary_statistics_id"], META["left_trait_id"], META["right_trait_id"], scale, display(covariance["estimate"]), display(covariance["standard_error"]), classification, artifact]
        )
    write_tsv(
        "heritability.tsv",
        ["relationship_id", "request_id", "method", "endpoint", "summary_statistics_id", "trait_id", "trait_type", "scale", "estimate", "standard_error", "classification", "native_artifact"],
        heritability_rows,
    )
    write_tsv(
        "genetic_covariance.tsv",
        ["relationship_id", "request_id", "method", "left_summary_statistics_id", "right_summary_statistics_id", "left_trait_id", "right_trait_id", "scale", "estimate", "standard_error", "classification", "native_artifact"],
        covariance_rows,
    )
    correlation = primary_correlation["correlation"]
    write_tsv(
        "genetic_correlation.tsv",
        ["relationship_id", "request_id", "method", "left_summary_statistics_id", "right_summary_statistics_id", "left_trait_id", "right_trait_id", "estimate", "standard_error", "z_score", "p_value", "classification", "native_artifact"],
        [[META["relationship_id"], META["request_id"], METHOD, META["left_summary_statistics_id"], META["right_summary_statistics_id"], META["left_trait_id"], META["right_trait_id"], display(correlation["estimate"]), display(correlation["standard_error"]), display(primary_correlation["z_score"]["estimate"]), display(primary_correlation["p_value"]["estimate"]), classification, "requests/{}/{}/native.observed.log".format(METHOD, META["request_id"])]],
    )
    native = dict(primary_diagnostics)
    native["estimates"] = native_by_scale

scales = list(native["estimates"].keys())
write_tsv(
    "diagnostics.tsv",
    ["relationship_id", "request_id", "metric", "value"],
    diagnostic_rows(classification, warning_codes, scales, native),
)

provenance = {
    "schema_version": "1.0",
    "request_id": META["request_id"],
    "primary_request_id": META["primary_request_id"],
    "request_name": META.get("request_name"),
    "method": METHOD,
    "relationship_id": META.get("relationship_id"),
    "endpoints": (
        {
            "left": {
                "summary_statistics_id": META["left_summary_statistics_id"],
                "trait_id": META["left_trait_id"],
                "trait_type": META["left_trait_type"],
                "population_prevalence": META.get("left_population_prevalence"),
                "sample_prevalence": META.get("left_sample_prevalence"),
            },
            "right": {
                "summary_statistics_id": META["right_summary_statistics_id"],
                "trait_id": META["right_trait_id"],
                "trait_type": META["right_trait_type"],
                "population_prevalence": META.get("right_population_prevalence"),
                "sample_prevalence": META.get("right_sample_prevalence"),
            },
        }
        if METHOD == "ldsc_rg"
        else {
            "unary": {
                "summary_statistics_id": META["summary_statistics_id"],
                "trait_id": META["trait_id"],
                "trait_type": META["trait_type"],
                "population_prevalence": META.get("population_prevalence"),
                "sample_prevalence": META.get("sample_prevalence"),
            }
        }
    ),
    "native_args": META["native_args"],
    "reference_bundle": META["reference_metadata"],
    "reference_validation": META["reference_validation"],
    "munging": {
        "keys": META["munging_keys"],
        "adapter_contract": "nfcore_gwas_canonical_v1_to_ldsc_sumstats_v1",
        "diagnostics": munging,
    },
    "classification": classification,
    "warnings": warning_codes,
    "native": native,
    "execution": {
        "normalizer_process": PROCESS_NAME,
        "normalizer_container": CONTAINER,
        "ldsc_container": LDSC_CONTAINER,
    },
    "tool": {"name": "ldsc", "native_version": native_version, "source_revision": "6c673952cee74bd5c57aef1555a03b1c015399a0"},
    "artifacts": {
        "native": ["native.{}.log".format(scale) for scale in scales],
        "normalized": ["heritability.tsv"] + (["genetic_covariance.tsv", "genetic_correlation.tsv"] if METHOD == "ldsc_rg" else []) + ["diagnostics.tsv", "provenance.json"],
    },
}
with open("provenance.json", "w", newline="") as handle:
    json.dump(provenance, handle, indent=2, sort_keys=True, allow_nan=False)
    handle.write(chr(10))

with open("versions.yml", "w", newline="") as handle:
    handle.write('"{}":'.format(PROCESS_NAME) + chr(10))
    handle.write("    python: {}".format(sys.version.split()[0]) + chr(10))
