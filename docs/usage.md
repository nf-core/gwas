# nf-core/gwas: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/gwas/usage](https://nf-co.re/gwas/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/gwas runs association, individual-level and summary-level heritability, and explicitly declared pairwise genetic-correlation analyses. A cohort manifest owns genotype facts; an analysis manifest links each individual-level trait analysis to one cohort; a summary-statistics manifest declares external or pipeline-generated summary results and their unary methods; and an optional relationship manifest binds ordered analysis or summary endpoints to pairwise methods.

> [!IMPORTANT]
> Genotypes must be prepared before you run the pipeline. The pipeline converts accepted genotype
> encodings into the formats required by its methods, but it does not perform genotype quality control.

## Relational manifest input

Every run supplies at least one input family: the linked `--cohort_manifest` and `--analysis_manifest`, `--summary_statistics_manifest`, or both. The linked cohort and analysis manifests must always be supplied together. `--relationship_manifest` is optional; when absent, no pair is inferred.

```bash
nextflow run nf-core/gwas \
    -r <VERSION> \
    -profile docker \
    --cohort_manifest ./cohorts.csv \
    --analysis_manifest ./analyses.csv \
    --outdir ./results
```

The cohort manifest owns reusable genotype facts. The analysis manifest owns one cohort–trait analysis and joins to exactly one cohort through `cohort_id`. Each `cohort_id` has one non-conflicting definition; each `analysis_id` is unique; and every analysis references a declared cohort. Cohort preparation and compatible matrices or predictions are reused, while outputs retain `analysis_id`.

### Cohort manifest fields

| Column         | Required | Description                                                                                                             |
| -------------- | -------- | ----------------------------------------------------------------------------------------------------------------------- |
| `cohort_id`    | Yes      | Unique, whitespace-free cohort identifier and analysis-manifest foreign key.                                            |
| `genome_build` | Yes      | `GRCh37` or `GRCh38`; selects build-specific harmonisation resources.                                                   |
| `ancestry`     | Yes      | Case-sensitive provenance label beginning with a letter or digit and containing only letters, digits, `_`, `.`, or `-`. |
| `pgen`         | By group | PLINK 2 genotype file; supply the complete `pgen`/`psam`/`pvar` group.                                                  |
| `psam`         | By group | PLINK 2 sample file.                                                                                                    |
| `pvar`         | By group | PLINK 2 variant file ending in `.pvar` or `.pvar.zst`.                                                                  |
| `bed`          | By group | PLINK 1 genotype file; supply the complete `bed`/`bim`/`fam` group.                                                     |
| `bim`          | By group | PLINK 1 variant file.                                                                                                   |
| `fam`          | By group | PLINK 1 sample file.                                                                                                    |
| `vcf`          | By group | One `.vcf`, `.vcf.gz`, or `.vcf.bgz` file; an index is not a manifest field.                                            |

Populate exactly one complete genotype representation on each row. PLINK 1 and VCF cohorts are converted once to canonical PLINK 2; PLINK 2 cohorts pass through.

### Analysis manifest fields

| Column                  | Required | Description                                                                                                                 |
| ----------------------- | -------- | --------------------------------------------------------------------------------------------------------------------------- |
| `analysis_id`           | Yes      | Unique, whitespace-free result identifier.                                                                                  |
| `cohort_id`             | Yes      | A `cohort_id` declared in the cohort manifest.                                                                              |
| `trait_id`              | Yes      | Whitespace-free trait identifier retained in provenance.                                                                    |
| `trait_type`            | Yes      | `quantitative` or `binary`; never inferred from values.                                                                     |
| `phenotype`             | Yes      | Existing headered phenotype file containing the selected trait.                                                             |
| `phenotype_column`      | Yes      | Header name of the selected trait column.                                                                                   |
| `control_value`         | Binary   | Source value recoded to `0`; required for binary traits, forbidden for quantitative traits, and distinct from `case_value`. |
| `case_value`            | Binary   | Source value recoded to `1`; required for binary traits and forbidden for quantitative traits.                              |
| `quant_covariates`      | No       | Existing headered quantitative-covariate file.                                                                              |
| `cat_covariates`        | No       | Existing headered categorical-covariate file.                                                                               |
| `association_methods`   | By row   | Optional comma-delimited selector: `plink2`, `regenie`, `gcta_fastgwa`, or `ldak_kvik`.                                     |
| `heritability_methods`  | By row   | Optional comma-delimited selector: `gcta_greml`, `gcta_greml_ldms`, `ldak_reml`, `ldak_he`, or `ldak_pcgc`.                 |
| `population_prevalence` | By route | Number strictly between `0` and `1`; valid only for binary heritability analyses and required by `ldak_pcgc`.               |
| `sample_prevalence`     | No       | Optional sample case fraction strictly between `0` and `1` for binary traits; forbidden for quantitative traits.          |

At least one unary method selector must be populated unless the analysis is referenced by a relationship row. Tokens are comma-delimited without spaces and may appear only once. Association-only, heritability-only and relationship-only analysis rows are valid.

### Summary-statistics manifest fields

The summary-statistics manifest has exactly sixteen columns and owns one stable `summary_statistics_id` per row. Each row declares exactly one mutually exclusive origin:

- An external result populates `source`, `source_mode` and `source_format`, leaves both producer fields blank, and declares its trait and provenance metadata.
- A pipeline-generated result populates `producer_analysis_id` and `producer_association_method`, leaves the three external source fields blank, and uses the exact deterministic ID `<producer_analysis_id>--<producer_association_method>`. Trait, build, ancestry, source method and prevalence are derived from the producer analysis and remain blank on the summary row.

| Column                        | Required        | Description                                                                                                                                       |
| ----------------------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| `summary_statistics_id`       | Yes             | Unique stable result identity. Internal results must use `<analysis_id>--<association_method>`.                                                    |
| `trait_id`                    | External        | Declared trait identity for an external result; derived for an internal result.                                                                   |
| `trait_type`                  | External        | `quantitative` or `binary` for an external result; derived for an internal result.                                                                |
| `source`                      | External        | Existing external `.tsv`, `.txt` or `.csv` table, optionally gzip-compressed.                                                                      |
| `source_mode`                 | External        | `raw` to harmonise through GWASLab or `canonical` to validate without reharmonising.                                                              |
| `source_format`               | External        | Explicit named GWASLab format for `raw`, or exactly `nfcore_gwas_canonical_v1` for `canonical`. Automatic format detection is not supported.       |
| `producer_analysis_id`        | Internal        | Declared analysis that produced the summary result.                                                                                               |
| `producer_association_method` | Internal        | Association method selected by that analysis.                                                                                                     |
| `genome_build`                | External        | `GRCh37` or `GRCh38`; derived from the producer cohort for an internal result.                                                                     |
| `ancestry`                    | External        | Researcher-declared provenance label; derived from the producer cohort for an internal result.                                                    |
| `source_method`               | External        | Program or method that produced the external table; derived from the producer association method for an internal result.                          |
| `source_release`              | No              | Optional external source or release token.                                                                                                        |
| `heritability_methods`        | By row          | Optional comma-delimited unary summary selector: `ldak_sumher` and/or `ldsc_h2`.                                                                  |
| `population_prevalence`       | Binary external | Optional population prevalence strictly between `0` and `1`; derived for an internal result.                                                      |
| `sample_prevalence`           | Binary external | Optional sample case fraction strictly between `0` and `1`; derived for an internal result.                                                       |
| `access_constraints`          | No              | Optional one-line, non-secret access or redistribution note retained in provenance; accepted for either origin.                                  |

An external canonical table must contain at least `SNPID`, `CHR`, `POS`, `EA`, `NEA`, `STATUS`, `EAF`, `BETA`, `SE`, `P` and `N`, with one non-empty variant row and a consistent tab-delimited field count. Raw and internal results converge on that same `nfcore_gwas_canonical_v1` contract before downstream summary methods consume them.

Each declared summary must either select a unary method or be referenced by a relationship. The primary unary request ID is deterministic:

```text
<method>--<summary_statistics_id>
```

### Relationship manifest fields

The optional relationship manifest has exactly eight columns. Endpoint slots are method-domain specific: `gcta_bivariate_reml` and `gcta_bivariate_reml_ldms` consume two analysis IDs, while `ldak_sumcors` and `ldsc_rg` consume two summary-statistics IDs. A row may select several methods for the same pair and may populate both endpoint domains only when each same-side summary is provably produced by the same-side analysis.

| Column                       | Required | Description                                                                                                                                          |
| ---------------------------- | -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| `relationship_id`            | Yes      | Unique, whitespace-free identifier for this exact populated endpoint binding.                                                                        |
| `left_analysis_id`           | GCTA     | Analysis ID for native trait 1.                                                                                                                       |
| `right_analysis_id`          | GCTA     | Analysis ID for native trait 2.                                                                                                                       |
| `left_summary_statistics_id` | Summary  | Summary-statistics ID bound to the ordered left endpoint for `ldak_sumcors` or `ldsc_rg`.                                                             |
| `right_summary_statistics_id` | Summary  | Summary-statistics ID bound to the ordered right endpoint for `ldak_sumcors` or `ldsc_rg`.                                                            |
| `relationship_methods`       | Yes      | Comma-delimited pairwise selectors: `gcta_bivariate_reml`, `gcta_bivariate_reml_ldms`, `ldak_sumcors`, and/or `ldsc_rg`.                              |
| `pair_quant_covariates`      | No       | Relationship-owned headered quantitative covariates beginning with `FID` and `IID`.                                                                  |
| `pair_cat_covariates`        | No       | Relationship-owned headered categorical covariates beginning with `FID` and `IID`.                                                                   |

The two populated endpoints must be different IDs and their declared trait IDs must differ. Different IDs may still represent related biological phenotypes; the pipeline does not police naming conventions. GCTA additionally requires two analysis IDs from one cohort. A reversed duplicate such as `height,disease` plus `disease,height` is invalid because it requests the same unordered binding twice. Left and right still matter and are recorded in every result and provenance file.

The pair phenotype is a deterministic full union of the two normalized endpoint sample sets in `FID`,`IID` order, with `NA` on a side where that trait is missing. Pair covariates belong to the relationship, not either endpoint analysis. The primary pair request ID is deterministic:

```text
gcta_bivariate_reml--<relationship_id>
gcta_bivariate_reml_ldms--<relationship_id>
```

Summary pair request IDs use the same rule:

```text
<method>--<relationship_id>
```

### Examples

[`assets/examples/relational/cohort_manifest.csv`](../assets/examples/relational/cohort_manifest.csv) is shared by the [minimal quantitative](../assets/examples/relational/analysis_manifest_quantitative.csv), [minimal binary](../assets/examples/relational/analysis_manifest_binary.csv), [association-only](../assets/examples/relational/analysis_manifest_association_only.csv), and [heritability-only](../assets/examples/relational/analysis_manifest_heritability_only.csv) examples. They use standard defaults and do not need `--method_options`.

```bash
nextflow run nf-core/gwas \
    -r <VERSION> \
    -profile docker \
    --cohort_manifest assets/examples/relational/cohort_manifest.csv \
    --analysis_manifest assets/examples/relational/analysis_manifest_quantitative.csv \
    --outdir results
```

The [heterogeneous multi-method manifest](../assets/examples/relational/analysis_manifest_heterogeneous.csv) pairs with [per-analysis GCTA and LDAK overrides](../assets/examples/relational/method_options_heterogeneous.json), including stageable resources:

```bash
nextflow run nf-core/gwas \
    -r <VERSION> \
    -profile docker \
    --cohort_manifest assets/examples/relational/cohort_manifest.csv \
    --analysis_manifest assets/examples/relational/analysis_manifest_heterogeneous.csv \
    --method_options assets/examples/relational/method_options_heterogeneous.json \
    --outdir results
```

The [relationship manifest example](../assets/examples/relational/relationship_manifest.csv) binds the heterogeneous quantitative and binary analyses and selects both GCTA pair estimators. For the deliberately small fixture, use the [namespaced method-options example](../assets/examples/relational/method_options_heterogeneous_bivariate.json), which retains the unary settings, uses one populated LDMS stratum, raises GCTA's native REML iteration allowance and deliberately drops the residual-covariance component. These are explicit example-specific choices, not pair-route defaults:

```bash
nextflow run nf-core/gwas \
    -r <VERSION> \
    -profile docker \
    --cohort_manifest assets/examples/relational/cohort_manifest.csv \
    --analysis_manifest assets/examples/relational/analysis_manifest_heterogeneous.csv \
    --relationship_manifest assets/examples/relational/relationship_manifest.csv \
    --method_options assets/examples/relational/method_options_heterogeneous_bivariate.json \
    --outdir results
```

The [mixed summary-statistics manifest](../assets/examples/relational/summary_statistics_manifest.csv) illustrates both origins: one internal PLINK 2 result from `heterogeneous_qt` and one already-canonical external result. The companion [summary relationship](../assets/examples/relational/relationship_manifest_summary.csv), [request options](../assets/examples/relational/method_options_summary.json), and [reference-catalog shape](../assets/examples/relational/reference_catalog.json) show the complete declaration surface. Replace every `/refs/...` value in the catalog with a locally available scientific reference before launching; the pipeline deliberately rejects unavailable paths and does not infer a bundle from ancestry.

### Advanced method options

`--method_options` is optional. The established form keeps its JSON root keyed by `analysis_id`; each value may contain `gcta`, `ldak` and/or `regenie`. It remains supported unchanged. A namespaced document places those same entries under `analyses`, summary unary settings under `unary_requests`, and relationship settings under `pair_requests`. Unlisted analysis settings receive their defaults. Every selected LDAK or LDSC summary request must explicitly choose a `reference_bundle_id`; nothing is inferred from the summary's ancestry label.

Both GCTA pair routes expose `native_args` as an array of individual non-file GCTA tokens on the deterministic request ID. A relationship's deterministic LDMS request additionally owns its matrix construction settings; it never inherits them from either endpoint's unary analysis:

```json
{
  "analyses": {},
  "pair_requests": {
    "gcta_bivariate_reml--height_disease": {
      "native_args": ["--reml-maxit", "500"]
    },
    "gcta_bivariate_reml_ldms--height_disease": {
      "ld_score_region_kb": 50,
      "ld_bins": 1,
      "ldms_maf_edges": [0, 0.5],
      "native_args": ["--reml-maxit", "500"]
    }
  }
}
```

The wrapper rejects whitespace or shell syntax, path separators, environment assignments, undeclared file-like values, file-bearing invocation mechanics such as `--keep` and `--extract`, wrapper-owned flags such as `--grm`, `--pheno`, `--out`, `--reml-bivar` and `--reml-bivar-prevalence`, and flags selecting another primary GCTA operation such as `--pca`. Arguments remain native scientific options: the pipeline records them and presents all resulting estimates; it does not choose a preferred result.

Summary requests use the same deterministic ownership boundary. LDAK receives exactly one staged `tagging_file`; LDSC receives separate staged `hapmap3_snplist`, `reference_ld_scores` and `regression_weights` roles. `native_args` may contain non-file scientific tokens only. Wrapper-owned operation, input, output and thread flags are rejected, as are any values that resemble undeclared files or paths.

For `ldak_sumher` and `ldak_sumcors`, the pipeline converts each distinct canonical summary once to LDAK's `Predictor A1 A2 Z n A1Freq` contract, with `A1` equal to the canonical effect allele and `Z = BETA / SE`; the canonical artifact and its provenance remain unchanged. Both routes use `--cutoff 0.01` unless a request explicitly supplies `--cutoff` or `--truncate`, and those two large-effect policies cannot be combined. SumCors initially accepts `LDAK-Thin`, `Uniform-GCTA` and `Human-Default` tagging bundles. Binary SumHer receives population prevalence and sample ascertainment only when both are declared. SumCors receives the two ordered prevalence/ascertainment pairs only when both endpoints are binary and all four values are present; mixed-trait and incomplete binary pairs run without liability arguments. LDAK's native ambiguous-variant exclusion and complete-summary checks remain enabled unless an accepted scientific override changes them.

LDSC H2 always retains an observed-scale invocation. A binary summary additionally receives a liability-scale invocation only when both `sample_prevalence` and `population_prevalence` are declared. LDSC RG follows the same rule for its ordered endpoints: the observed-scale result is always retained, and a second invocation supplies prevalence only when at least one endpoint is binary and every binary endpoint declares both values. Quantitative endpoints in that mixed invocation use LDSC's native `nan` prevalence placeholder. Request-level `native_args` cannot override either prevalence flag.

```json
{
  "unary_requests": {
    "ldsc_h2--height--plink2": {
      "reference_bundle_id": "ldsc_eur"
    }
  },
  "pair_requests": {
    "ldak_sumcors--height_disease": {
      "reference_bundle_id": "ldak_thin_eur"
    },
    "ldsc_rg--height_disease": {
      "reference_bundle_id": "ldsc_eur"
    }
  }
}
```

Each selected entity–method binding creates one deterministic primary request. Additional configurations are complete, independent named requests; they do not inherit omitted settings from the primary. Declare `primary_request_id` and `request_name`, and key the object by `<primary_request_id>--<request_name>`:

```json
{
  "unary_requests": {
    "ldsc_h2--height--plink2": {
      "reference_bundle_id": "ldsc_eur"
    },
    "ldsc_h2--height--plink2--alternate": {
      "primary_request_id": "ldsc_h2--height--plink2",
      "request_name": "alternate",
      "reference_bundle_id": "ldsc_eur_alternate"
    }
  }
}
```

### Reference catalog

`--reference_catalog` is a JSON object with optional `ldsc` and `ldak` families. Bundle IDs are globally unique across both families.

```json
{
  "ldsc": {
    "ldsc_eur": {
      "genome_build": "GRCh37",
      "ancestry": "EUR",
      "variant_id_system": "rsid",
      "hapmap3_snplist": "/refs/w_hm3.snplist",
      "reference_ld_scores": "/refs/eur_w_ld_chr/",
      "regression_weights": "/refs/eur_w_ld_chr/"
    }
  },
  "ldak": {
    "ldak_thin_eur": {
      "genome_build": "GRCh37",
      "ancestry": "EUR",
      "variant_id_system": "rsid",
      "model": "LDAK-Thin",
      "tagging_file": "/refs/ldak-thin.tagging"
    }
  }
}
```

The catalog may also declare a SHA-256 digest beside each role. Preflight checks the document shape, family, required roles, digest syntax and path availability. In the first release, `genome_build`, `ancestry`, `variant_id_system` and `model` are recorded provenance rather than a pipeline certification of scientific compatibility. The user owns reference selection.

| GCTA option          | Type and default                      | Consumer and constraints                                                                                |
| -------------------- | ------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| `grm_maf`            | Number or `null`; `null`              | Dense `gcta_greml` matrix filter; `0` to `0.5` inclusive when set.                                      |
| `grm_extract`        | Resource path or absent; absent       | Predictor list staged for dense `gcta_greml` matrix construction.                                       |
| `reml_no_constrain`  | Boolean; `false`                      | `gcta_greml` and `gcta_greml_ldms`; disables variance-component constraints.                            |
| `sparse_cutoff`      | Number; `0.05`                        | `gcta_fastgwa` only; `0` to `1` inclusive.                                                              |
| `ld_score_region_kb` | Positive integer; `200`               | `gcta_greml_ldms` only; LD-score window in kilobases.                                                   |
| `ld_bins`            | Positive integer; `4`                 | `gcta_greml_ldms` only; number of LD-score strata.                                                      |
| `ldms_maf_edges`     | Number array; `[0,0.01,0.05,0.2,0.5]` | `gcta_greml_ldms` only; strictly increasing, at least two values, beginning at `0` and ending at `0.5`. |

The same three LDMS setting names are accepted under a `gcta_bivariate_reml_ldms--<relationship_id>` pair request. Their defaults are `200`, `4` and `[0,0.01,0.05,0.2,0.5]`. The resolved values become part of the matrix reuse key, so a unary and pair request reuse one MGRM family only when cohort, genotype input and all three settings agree exactly.

| LDAK option          | Type and default                | Consumer and constraints                                                                                                      |
| -------------------- | ------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| `model`              | String; `human_default`         | Kinship routes; `human_default` or `custom`.                                                                                  |
| `power`              | Number; `-0.25`                 | Kinship routes; `-2` to `0`. `human_default` fixes `-0.25`; select `custom` for another value.                                |
| `weights_policy`     | String; `equal`                 | Kinship routes; `equal` passes `--ignore-weights YES`, `default` retains native weighting, and `provided` requires `weights`. |
| `weights`            | Resource path or absent; absent | Kinship routes (`ldak_reml`, `ldak_he`, `ldak_pcgc`); accepted exactly with `weights_policy: provided`.                       |
| `relatedness_filter` | Boolean; `false`                | Kinship routes; optionally derives an unrelated subset.                                                                       |
| `kvik_step1_subset`  | String; `all`                   | `ldak_kvik` only; `all`, `thin_common`, or `provided`.                                                                        |
| `predictor_extract`  | Resource path or absent; absent | `ldak_kvik` only; required exactly with `kvik_step1_subset: provided`.                                                        |

| REGENIE option      | Type and default         | Consumer and constraints                                                                             |
| ------------------- | ------------------------ | ---------------------------------------------------------------------------------------------------- |
| `step1_bsize`       | Positive integer; `1000` | Step 1 fitted-model block size; participates in prediction-reuse identity.                           |
| `firth`             | Boolean; `true`          | Binary traits only; enable Firth fallback in Step 2.                                                 |
| `firth_approx`      | Boolean; `true`          | Binary traits only; requires `firth` when explicitly enabled.                                        |
| `firth_p_threshold` | Number; `0.01`           | Binary traits only; greater than `0` and at most `1`, and requires `firth` when explicitly supplied. |
| `min_mac`           | Number or `null`; `null` | Optional Step 2 minimum minor allele count; `null` leaves REGENIE's built-in behavior in effect.     |

For example, this changes the fitted-model block size and Step 2 policy for one binary REGENIE analysis while every unlisted analysis retains the defaults:

```json
{
  "disease": {
    "regenie": {
      "step1_bsize": 2000,
      "firth": false,
      "min_mac": 10
    }
  }
}
```

Resource paths are staged and their contents participate in matrix or prediction reuse identity. `gcta_grm_parts` is operational partitioning and remains configuration/profile-only, never a method option.

The authoritative structural contracts are [`schema_cohort_manifest.json`](../assets/schema_cohort_manifest.json), [`schema_analysis_manifest.json`](../assets/schema_analysis_manifest.json), [`schema_summary_statistics_manifest.json`](../assets/schema_summary_statistics_manifest.json), [`schema_relationship_manifest.json`](../assets/schema_relationship_manifest.json), and [`schema_reference_catalog.json`](../assets/schema_reference_catalog.json); origin-, relationship-, request-, method-, trait- and resource-aware diagnostics come from the central preflight validator.

### Validation diagnostics

Structural failures name the manifest and invalid column. Cross-row preflight failures additionally name the CSV row and `cohort_id`, `analysis_id`, `summary_statistics_id` or `relationship_id`: examples include an incomplete genotype group, conflicting duplicates, an invalid or mixed summary origin, an internal producer mismatch, orphan endpoints, same-endpoint or same-trait pairs, cross-cohort GCTA pairs, reversed duplicates, unknown or repeated method tokens, invalid binary coding, and route-inapplicable prevalence. Method-options and reference-catalog failures name the document, request or bundle ID, qualified option or resource role, and reason before task submission.

### Phenotype normalisation

The pipeline normalises each selected phenotype to a common two-identifier-plus-trait layout. Binary
source values matching `control_value` and `case_value` become `0` and `1`; missing values and unmatched
binary values become `NA`. Quantitative values must be numeric. Quantitative and categorical covariates
remain distinct for tools with separate native interfaces.

### Run-level defaults

| Parameter or behaviour            | Default and rationale                                                                                                                                                                                                                                                                      |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| PLINK 2 binary association        | Firth fallback is enabled so separated or sparse binary-trait tests can still produce estimates. Binary phenotypes are passed with `--1` because the common normalised coding is `0`/`1`/`NA`; covariates are variance-standardised to prevent numerical failure when their scales differ. |
| `--regenie_step1_mode`            | `standard`, the simplest one-task Step 1. Use `chunked` with `--regenie_step1_jobs` when a large cohort needs REGENIE's split-L0/run-L0/run-L1 execution family.                                                                                                                           |
| `--regenie_lowmem`                | `true`, keeping Step 1's temporary prediction blocks in the task work directory to reduce memory use.                                                                                                                                                                                      |
| REGENIE scientific method options | Per-analysis `regenie.*` defaults enable approximate Firth fallback below `0.01` for binary traits and leave `min_mac` unset so REGENIE's own versioned policy applies.                                                                                                                    |
| GWASLab reference parameters      | Unset. Every association output is still standardised; reference-dependent allele checks, rsID assignment and strand inference run only when you provide the corresponding build-specific FASTA or VCF resource.                                                                           |
| Save controls                     | Off. Intermediates stay out of the results directory unless explicitly requested, avoiding unexpectedly large published output.                                                                                                                                                            |

The four opt-in save controls are:

- `--save_prepared_genotypes`: publish PLINK 2 bundles that the pipeline converted under `genotypes/<cohort_id>/`.
- `--save_normalised_phenotypes`: publish headered normalised phenotype and covariate files under `phenotypes/<analysis_id>/`.
- `--save_relatedness_matrices`: publish merged GCTA or LDAK matrix bundles under `quality_control/relatedness_matrices/<key>/`.
- `--save_association_predictions`: publish reusable REGENIE and LDAK-KVIK Step 1 bundles under `intermediates/association_predictions/<method>/<analysis_id>/`.

See the [output documentation](output.md) for the exact files and publication exceptions.

## Running the pipeline

The typical command for running the pipeline is:

```bash
nextflow run nf-core/gwas \
    -r <VERSION> \
    -profile docker \
    --cohort_manifest ./cohorts.csv \
    --analysis_manifest ./analyses.csv \
    --outdir ./results
```

`-r <VERSION>` is strongly recommended for reproducibility. Choose the software profile appropriate for your system; see [`-profile`](#-profile).

The pipeline creates the following in your launch directory:

```bash
work                # Nextflow task work directories
<OUTDIR>            # Published results selected by --outdir
.nextflow_log       # Nextflow execution log
```

For repeated runs, place pipeline parameters in a YAML or JSON file and supply it with `-params-file`:

```bash
nextflow run nf-core/gwas -r <VERSION> -profile docker -params-file params.yaml
```

```yaml title="params.yaml"
cohort_manifest: "./cohorts.csv"
analysis_manifest: "./analyses.csv"
relationship_manifest: "./relationships.csv"
outdir: "./results/"
```

You can also generate parameter files with [nf-core/launch](https://nf-co.re/launch).

> [!WARNING]
> Do not use `-c <file>` to specify pipeline parameters. Custom configuration files loaded with `-c` are for process resources, infrastructure settings and module arguments; use CLI flags or `-params-file` for parameters.

## Updating

Nextflow caches pipeline code after the first run. Update the cached copy with:

```bash
nextflow pull nf-core/gwas
```

## Reproducibility

Pin the pipeline revision with `-r`, retain the exact supplied manifests, reference catalog, method-options document and parameter file, and archive `pipeline_info/` with your results. Reusing the same revision, inputs and parameters also lets `-resume` recover cached tasks.

Find published version tags on the [nf-core/gwas versions page](https://github.com/nf-core/gwas/tags). The run's pipeline and tool versions are recorded in the published reports described in [Pipeline information](output.md#pipeline-information).

## Core Nextflow arguments

> [!NOTE]
> Nextflow options use one hyphen; pipeline parameters use two.

### `-profile`

Use `-profile` to select configuration for your software or compute environment. Multiple profiles are comma-separated and applied from left to right, for example `-profile test,docker`.

The generic software profiles include `docker`, `singularity`, `apptainer`, `podman`, `shifter`, `charliecloud`, `conda`, `mamba` and `wave`. Docker or Singularity/Apptainer is recommended for reproducibility. If no profile is supplied, every tool must already be available on `PATH`.

The pipeline also loads institutional profiles from [nf-core/configs](https://github.com/nf-core/configs#documentation).

### `-resume`

Add `-resume` when restarting a run. Nextflow reuses cached tasks whose inputs, code and configuration have not changed. You can resume a named run with `-resume <run-name>`; use `nextflow log` to list run names.

### `-c`

Use `-c <file>` to load a Nextflow configuration file for infrastructure, executors and process resources. See the [nf-core configuration documentation](https://nf-co.re/docs/running/configuration).

## Custom configuration

### Resource requests

The pipeline assigns default CPU, memory and time requests by process label and retries selected resource-related failures with larger requests. Use [`--max_cpus`, `--max_memory` and `--max_time`](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) to cap requests, or use a custom configuration file to [customise individual processes](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources).

### Custom containers

To override a tool container or Conda environment, follow the nf-core guidance for [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions). Record overrides because they change the software provenance of your run.

### Custom tool arguments

Use a custom configuration file to add supported process-specific arguments as described in [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments). Do not use this mechanism to reach methods outside the documented selector lists.

### nf-core/configs

If a configuration is useful across your institution, test it locally and consider contributing it to [nf-core/configs](https://github.com/nf-core/configs). Include the corresponding documentation and profile registration in that repository.

## Running in the background

Nextflow must keep running to submit and supervise jobs. Use `-bg`, a detached `screen` or `tmux` session, or a cluster job appropriate to your infrastructure. With `-bg`, Nextflow writes its console output to a log file.

## Nextflow memory requirements

The Nextflow Java process can require substantial memory when coordinating a large run. Limit its heap with `NXF_OPTS`, for example:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```

## Troubleshooting

### Start with the reported phase

| What you see                                             | Meaning                                                                                              | First evidence to inspect                                                                                                |
| -------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| An error naming a manifest, CSV row, identity and field            | Preflight rejected input before task submission                                                      | Read the complete grouped error report, fix every listed row and relaunch.                                               |
| `Method-options document ... request_id ... option ...`            | A deterministic or named request is incomplete, invalid or conflicts with the invocation firewall    | Read the request namespace, request ID, option and reason; then correct the complete request configuration.              |
| `Reference catalog ... reference_bundle_id ... field ...`          | A bundle is malformed, duplicated across families or has an unavailable required resource            | Correct the named bundle role; compatibility metadata is recorded but not scientifically certified by the pipeline.    |
| `Process ... terminated with an error`                             | A task was submitted and failed                                                                      | Inspect `.nextflow.log` and the task's `.command.err`, `.command.out` and `.command.log` in the reported work directory. |

### Manifest validation failed

Read the reported CSV path, row number, entity ID, field name and reason from left to right. The validator reports all linked-manifest errors it can find in one launch, so correct every bullet before rerunning.

| Diagnostic fragment                                                                                                       | Supported action                                                                                                                                                                                            |
| ------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `header row 1 does not match the mandatory ...-column input contract`                                                     | Start from a shipped example and restore missing columns; remove unexpected or repeated column names. Optional values may be blank, but their columns must remain present. Column order is not significant. |
| `no genotype group is populated`, `a second genotype group is populated` or `genotype group ... is only partly populated` | Populate exactly one complete representation: `pgen`/`psam`/`pvar`, `bed`/`bim`/`fam` or `vcf`.                                                                                                             |
| `duplicate cohort_id` or `duplicate analysis_id`                                                                          | Keep one row for each identity. Conflicting duplicates also report the fields that differ.                                                                                                                  |
| `undefined cohort_id`                                                                                                     | Make the analysis row's `cohort_id` match one cohort-manifest identity exactly.                                                                                                                             |
| `unknown method`, `listed more than once` or `row selects no method`                                                      | Use each documented selector at most once and populate at least one of `association_methods` or `heritability_methods`. See [Analysis manifest fields](#analysis-manifest-fields).                          |
| A binary `case_value` or `control_value` diagnostic                                                                       | Supply both distinct source codes for a binary trait. Remove both from a quantitative row. The pipeline does not infer `trait_type` from phenotype values.                                                  |
| A `population_prevalence` diagnostic                                                                                      | Use it only for a binary heritability analysis whose selected estimator consumes it, and provide it for `ldak_pcgc`. See [Analysis manifest fields](#analysis-manifest-fields).                             |
| A summary origin or deterministic identity diagnostic                                                                     | Populate exactly one complete external or producer origin. Internal IDs must be `<producer_analysis_id>--<producer_association_method>`. See [Summary-statistics manifest fields](#summary-statistics-manifest-fields). |
| An undefined or self-paired summary endpoint                                                                               | Declare each referenced summary ID, use distinct endpoint and trait IDs on the two sides, and ensure any combined analysis/summary side has recorded producer correspondence.                              |

### Method-options validation failed

The diagnostic names the document, entity or request ID, fully qualified option and reason. Fix malformed JSON; use the `gcta`, `ldak` and `regenie` analysis families; and use only the `analyses`, `unary_requests` and `pair_requests` namespaces described under [Advanced method options](#advanced-method-options). A request may configure only a method already selected by its summary or relationship declaration. Named additions must provide their own complete settings, including a reference bundle for LDAK or LDSC.

Resolve resource paths from the launch environment. `weights_policy: provided` requires `weights`, and `kvik_step1_subset: provided` requires `predictor_extract`; resources are also rejected when supplied under an incompatible policy. Operational settings such as process resources and tool threads belong in run or profile configuration, not `--method_options`.

### A process failed

1. Copy the process name and work-directory hash from the terminal or `.nextflow.log`.
2. Inspect `.command.err`, `.command.out` and `.command.log` in that task directory; `.command.sh` records the executed command.
3. Correct the reported input, configuration, container, filesystem or resource problem.
4. Relaunch the same command with `-resume`; do not delete `work/` or `.nextflow/` first.

Selected resource failures are retried with larger requests, subject to `--max_cpus`, `--max_memory` and `--max_time`. See [Resource requests](#resource-requests) and the general [nf-core troubleshooting guide](https://nf-co.re/docs/running/troubleshooting).

### Expected results are missing

- Confirm that the analysis, summary-statistics or relationship row selected the method whose result you expected.
- Check the route-specific paths in the [output documentation](output.md).
- Prepared genotypes, normalised phenotypes and covariates, relatedness matrices and REGENIE predictions are unpublished by default; enable the corresponding save control before expecting those directories.
- GWASLab reference parameters are optional. Standardised association output is still produced without them, but reference-dependent allele checks, flips, rsID assignment and strand inference are not. `genome_build`, not `ancestry`, selects the build-specific resources.
- Published intermediates are retention outputs, not importable cross-run caches. Supported reuse requires retained Nextflow work and `-resume`.

### `-resume` did not reuse work

Nextflow reuses a task only when its inputs, pipeline code and relevant configuration still match a retained cache entry. Keep the original `work/` directory and `.nextflow/` cache, rerun from the same launch context where practical, and use `nextflow log` to find a run name for `-resume <run-name>`. A changed input file, manifest, method option, pipeline revision, process configuration or missing work directory can require recomputation. See [`-resume`](#-resume).

### Get help

Before requesting help, retain the pipeline version, exact launch command with secrets removed, `.nextflow.log`, the failing process name and relevant `.command.*` files, the reported manifest row or a minimal redacted reproducer, executor/profile and container runtime, and whether the original work directory remains available for `-resume`.

Use the [nf-core troubleshooting guide](https://nf-co.re/docs/running/troubleshooting), open a [GitHub issue](https://github.com/nf-core/gwas/issues) with a reproducible report, or ask in the nf-core [`#gwas` Slack channel](https://nfcore.slack.com/channels/gwas).
