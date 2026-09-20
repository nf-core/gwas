# nf-core/gwas: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/gwas/usage](https://nf-co.re/gwas/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/gwas runs association, individual-level and summary-level heritability, and explicitly declared pairwise genetic-correlation analyses. A cohort manifest owns genotype facts; an analysis manifest links each individual-level trait analysis to one cohort; a summary-statistics manifest declares external, pipeline-generated or meta-analysis-derived summary results and their unary methods; and an optional relationship manifest binds ordered analysis or summary endpoints to pairwise methods.

> [!IMPORTANT]
> Genotypes must be prepared before you run the pipeline. The pipeline converts accepted genotype
> encodings into the formats required by its methods, but it does not perform genotype quality control.

## Relational manifest input

Every run supplies at least one input family: the linked `--cohort_manifest` and `--analysis_manifest`, `--summary_statistics_manifest`, or both. The linked cohort and analysis manifests must always be supplied together. `--relationship_manifest` is optional; when absent, no pair is inferred.

The JSON schemas own column names, required fields, types, enumerations and file existence. Schema-optional columns may be omitted entirely or supplied with blank cells. Unexpected and repeated column names are rejected; quoted CSV headers and values are supported, and column order is not significant.

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

| Column             | Required | Description                                                                                                                                                                                                                              |
| ------------------ | -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `cohort_id`        | Yes      | Unique, whitespace-free cohort identifier and analysis-manifest foreign key.                                                                                                                                                             |
| `genome_build`     | Yes      | `GRCh37` or `GRCh38`; selects build-specific harmonisation resources.                                                                                                                                                                    |
| `ancestry`         | Yes      | Case-sensitive provenance label beginning with a letter or digit and containing only letters, digits, `_`, `.`, or `-`.                                                                                                                  |
| `pgen`             | By group | PLINK 2 genotype file; supply the complete `pgen`/`psam`/`pvar` group.                                                                                                                                                                   |
| `psam`             | By group | PLINK 2 sample file.                                                                                                                                                                                                                     |
| `pvar`             | By group | PLINK 2 variant file ending in `.pvar`. A Zstandard-compressed `.pvar.zst` is not accepted: it requires PLINK 2's `vzs` modifier, which the PLINK 2, GCTA and REGENIE consumers here do not pass.                                        |
| `bed`              | By group | PLINK 1 genotype file; supply the complete `bed`/`bim`/`fam` group.                                                                                                                                                                      |
| `bim`              | By group | PLINK 1 variant file.                                                                                                                                                                                                                    |
| `fam`              | By group | PLINK 1 sample file.                                                                                                                                                                                                                     |
| `genotype_view_id` | No       | Immutable identity of the supplied genotype view: 16–64 lowercase hexadecimal characters, optionally prefixed by a backend namespace such as `biofuse:`. When absent, the pipeline digests the supplied bytes once per cohort in a task. |

Populate exactly one complete genotype representation on each row, and give the members of a multi-file group one shared basename stem — every PLINK, GCTA and REGENIE consumer addresses a fileset by a single prefix.

The supplied representation is preserved. A PLINK 1 or PLINK 2 cohort is used exactly as given and nothing is converted for it. A PLINK 1 BED/BIM/FAM view is derived from a PGEN cohort once, and only when a selected method reads PLINK 1 — the LDAK estimators, LDAK-KVIK, the LD- and MAF-stratified GCTA routes and the MPH routes. That projection takes hard calls at an explicit `--hard-call-threshold 0.1` and drops dosage and phase; every cohort's `genotypes/<cohort_id>/<cohort_id>.genotype_view.json` records which of these applied to it.

### Analysis manifest fields

| Column                  | Required | Description                                                                                                                                                                |
| ----------------------- | -------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `analysis_id`           | Yes      | Unique, whitespace-free result identifier.                                                                                                                                 |
| `cohort_id`             | Yes      | A `cohort_id` declared in the cohort manifest.                                                                                                                             |
| `trait_id`              | Yes      | Whitespace-free trait identifier retained in provenance.                                                                                                                   |
| `trait_type`            | Yes      | `quantitative` or `binary`; never inferred from values.                                                                                                                    |
| `phenotype`             | Yes      | Existing headered phenotype file containing the selected trait.                                                                                                            |
| `phenotype_column`      | Yes      | Header name of the selected trait column.                                                                                                                                  |
| `control_value`         | Binary   | Source value recoded to `0`; required for binary traits, forbidden for quantitative traits, and distinct from `case_value`.                                                |
| `case_value`            | Binary   | Source value recoded to `1`; required for binary traits and forbidden for quantitative traits.                                                                             |
| `quant_covariates`      | No       | Existing headered quantitative-covariate file.                                                                                                                             |
| `cat_covariates`        | No       | Existing headered categorical-covariate file.                                                                                                                              |
| `association_methods`   | By row   | Optional comma-delimited selector: `regenie`, `gcta_fastgwa`, or `ldak_kvik`.                                                                                              |
| `heritability_methods`  | By row   | Optional comma-delimited selector: `gcta_greml`, `gcta_greml_ldms`, `ldak_reml`, `ldak_he`, `ldak_pcgc`, `ldak_fast_he`, `ldak_fast_pcgc`, `mph_reml`, or `mph_reml_ldms`. |
| `population_prevalence` | By route | Number strictly between `0` and `1`; valid only for binary heritability analyses and required by `ldak_pcgc` and `ldak_fast_pcgc`.                                         |
| `sample_prevalence`     | No       | Optional sample case fraction strictly between `0` and `1` for binary traits; forbidden for quantitative traits.                                                           |

At least one unary method selector must be populated unless the analysis is referenced by a relationship row. Tokens are comma-delimited without spaces and may appear only once. Association-only, heritability-only and relationship-only analysis rows are valid.

### Summary-statistics manifest fields

The summary-statistics manifest owns one stable `summary_statistics_id` per row. Each row declares exactly one mutually exclusive origin:

- An external result populates `source` and `source_format`, leaves both producer fields blank, and declares its trait metadata. Every external source passes through GWASLab; use the `gwaslab` format for a pre-harmonised GWASLab table.
- A pipeline-generated result populates `producer_analysis_id` and `producer_association_method`, leaves the two external source fields blank, and uses the exact deterministic ID `<producer_analysis_id>--<producer_association_method>`. Trait, build, ancestry, source method and prevalence are derived from the producer analysis and remain blank on the summary row.

| Column                        | Required        | Description                                                                                                                                                                                                                                                                                                                                 |
| ----------------------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `summary_statistics_id`       | Yes             | Unique stable result identity. Internal results must use `<analysis_id>--<association_method>`.                                                                                                                                                                                                                                             |
| `trait_id`                    | External        | Declared trait identity for an external result; derived for an internal result.                                                                                                                                                                                                                                                             |
| `trait_type`                  | External        | `quantitative` or `binary` for an external result; derived for an internal result.                                                                                                                                                                                                                                                          |
| `source`                      | External        | Existing external `.tsv`, `.txt` or `.csv` table, optionally gzip-compressed.                                                                                                                                                                                                                                                               |
| `source_format`               | External        | Named GWASLab input format; use `gwaslab` for a pre-harmonised GWASLab table. GWASLab's `auto*` formats apply the effect-allele orientation their name declares (for example `auto_2` treats `A2` as the effect allele) without verifying it against the file; the assumption GWASLab applied is recorded in the published `*.gwaslab.log`. |
| `producer_analysis_id`        | Internal        | Declared analysis that produced the summary result.                                                                                                                                                                                                                                                                                         |
| `producer_association_method` | Internal        | Association method selected by that analysis.                                                                                                                                                                                                                                                                                               |
| `genome_build`                | External        | `GRCh37` or `GRCh38`; derived from the producer cohort for an internal result.                                                                                                                                                                                                                                                              |
| `ancestry`                    | External        | Researcher-declared provenance label; derived from the producer cohort for an internal result.                                                                                                                                                                                                                                              |
| `source_method`               | External        | Program or method that produced the external table; derived from the producer association method for an internal result.                                                                                                                                                                                                                    |
| `source_release`              | No              | Optional external source or release token.                                                                                                                                                                                                                                                                                                  |
| `heritability_methods`        | By row          | Optional comma-delimited unary summary selector: `ldak_sumher` and/or `ldsc_h2`.                                                                                                                                                                                                                                                            |
| `population_prevalence`       | Binary external | Optional population prevalence strictly between `0` and `1`; derived for an internal result.                                                                                                                                                                                                                                                |
| `sample_prevalence`           | Binary external | Optional sample case fraction strictly between `0` and `1`; derived for an internal result.                                                                                                                                                                                                                                                 |
| `access_constraints`          | No              | Optional one-line, non-secret access or redistribution note associated with the declared summary result; accepted for either origin.                                                                                                                                                                                                        |

Every external and internal source crosses `GWASLAB_HARMONIZE`. The process emits the pipeline-standard GWASLab table directly; the pipeline does not add a second serializer, schema validator or provenance sidecar after GWASLab.

Each declared base summary must select a unary method, be referenced by a relationship, or be a parent of a meta-analysis. A meta-analysis row itself declares work and can also select downstream unary or pair methods. The primary unary request ID is deterministic:

```text
<method>--<summary_statistics_id>
```

The third origin uses `source_summary_statistics_ids`, `meta_analysis_models` and `min_studies`; see [Common-variant meta-analysis](#common-variant-meta-analysis) for ordered parent selection and inherited metadata.

### Relationship manifest fields

The optional relationship manifest uses method-domain-specific endpoint slots: `gcta_bivariate_reml`, `gcta_bivariate_reml_ldms`, `gcta_bivariate_he`, `gcta_bivariate_he_ldms`, `mph_bivariate_reml` and `mph_bivariate_reml_ldms` consume two analysis IDs, while `ldak_sumcors` and `ldsc_rg` consume two summary-statistics IDs. A row may select several methods for the same pair and may populate both endpoint domains only when each same-side summary is provably produced by the same-side analysis.

| Column                        | Required | Description                                                                                                                                                                                                              |
| ----------------------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `relationship_id`             | Yes      | Unique, whitespace-free identifier for this exact populated endpoint binding.                                                                                                                                            |
| `left_analysis_id`            | Analysis | Analysis ID for native trait 1.                                                                                                                                                                                          |
| `right_analysis_id`           | Analysis | Analysis ID for native trait 2.                                                                                                                                                                                          |
| `left_summary_statistics_id`  | Summary  | Summary-statistics ID bound to the ordered left endpoint for `ldak_sumcors` or `ldsc_rg`.                                                                                                                                |
| `right_summary_statistics_id` | Summary  | Summary-statistics ID bound to the ordered right endpoint for `ldak_sumcors` or `ldsc_rg`.                                                                                                                               |
| `relationship_methods`        | Yes      | Comma-delimited pairwise selectors: `gcta_bivariate_reml`, `gcta_bivariate_reml_ldms`, `gcta_bivariate_he`, `gcta_bivariate_he_ldms`, `mph_bivariate_reml`, `mph_bivariate_reml_ldms`, `ldak_sumcors`, and/or `ldsc_rg`. |
| `pair_quant_covariates`       | No       | Relationship-owned headered quantitative covariates beginning with `FID` and `IID`.                                                                                                                                      |
| `pair_cat_covariates`         | No       | Relationship-owned headered categorical covariates beginning with `FID` and `IID`.                                                                                                                                       |

The two populated endpoints must be different IDs and their declared trait IDs must differ. Different IDs may still represent related biological phenotypes; the pipeline does not police naming conventions. Every analysis-endpoint method additionally requires two analysis IDs from one cohort. A reversed duplicate such as `height,disease` plus `disease,height` is invalid because it requests the same unordered binding twice. Left and right still matter and are retained in request identity and native execution order.

`gcta_bivariate_he` and `gcta_bivariate_he_ldms` add GCTA's `--HEreg-bivar` Haseman-Elston cross-product estimator alongside `gcta_bivariate_reml` and `gcta_bivariate_reml_ldms`, using one dense GRM or the same LD- and MAF-stratified `--mgrm` family respectively. Both are explicit-GRM moment estimators: HE-CP only makes the _fitting_ stage cheaper than REML, and it still requires the same full dense (or LDMS-stratified) matrix construction, storage and I/O. Document and treat them as a deterministic moment reference and sensitivity analysis, not as a matrix-free or more scalable route.

The two HE selectors accept quantitative pairs only; a binary or mixed-endpoint pair is rejected before execution. `gcta_bivariate_reml` and `gcta_bivariate_reml_ldms` remain the supported route for binary and mixed pairs because they carry an explicit prevalence and liability-scale contract — they are not a slow fallback for those traits. Neither HE selector supports covariate adjustment: GCTA 1.94.1 lists `--qcovar` and `--covar` under "Accepted options" for `--HEreg-bivar` and then silently ignores them, so a relationship that declares `pair_quant_covariates` or `pair_cat_covariates` together with an HE selector is a validation error before execution; use `gcta_bivariate_reml` or `gcta_bivariate_reml_ldms` when covariate adjustment is required.

HE-CP regresses the cross-trait product on the lower triangle of the GRM only, with the left trait on the row member and the right trait on the column member, so the declared left/right orientation changes the point estimate — this is why a reversed duplicate relationship remains a validation error for the HE routes as well as the REML routes. Selecting `gcta_bivariate_he` (or `gcta_bivariate_he_ldms`) alongside `gcta_bivariate_reml` (or `gcta_bivariate_reml_ldms`) for the same relationship builds the required GRM or MGRM family only once: the matrix reuse key is derived from the cohort, the genotype bundle and the declared matrix settings, never from the method token.

`mph_bivariate_reml` and `mph_bivariate_reml_ldms` fit the same declared, oriented pair by MPH's REML, over one MPH relationship matrix or over the LD- and MAF-stratified MPH matrix family respectively. They are relationship methods in exactly the sense the GCTA pair routes are: they take two analysis endpoints from one cohort, their request identity is `<method>--<relationship_id>`, and their results publish under `requests/<method>/<request_id>/`. Both accept quantitative pairs only. A binary or mixed pair is refused at ingress, and the diagnostic names `gcta_bivariate_reml` and `gcta_bivariate_reml_ldms`, which carry the prevalence and liability-scale contract; MPH has no prevalence option at all.

The one thing to settle before selecting an MPH pair route beside a GCTA one is that they do not fit the same individuals. MPH fits both traits only on the individuals who have both traits and every named covariate observed. GCTA's bivariate REML does not: given the same partially overlapping pair it keeps every individual and fits an unbalanced design. Measured on a 200-individual cohort with the right trait blanked for 80 of them, MPH reported `Non-missing analysis set contains 120 individuals.` while GCTA reported `200 non-missing phenotypes for trait #1 and 120 for trait #2` and fitted 320 observations. An MPH pair result and a GCTA pair result over a partially overlapping relationship are therefore not matched models, and the difference between the two estimates is not attributable to the estimator alone. The native log and solver trace report the analysis-set size.

The declared left/right order is written into MPH's `--trait_names` in that order and is the execution orientation. MPH then labels the pair in its own output in the reverse order — `--trait_names A,B` produces rows labelled `trait_x = B`, `trait_y = A` — consistently, in `native.mq.cor.csv` and in the cross block of `native.mq.vc.csv`. The pipeline matches those rows by trait name and never by position, and records what it saw as `correlations.native_trait_pair_labels` and `correlations.native_pair_order`. If you read the native files yourself, read them by trait name: the column order is not the declared order.

Selecting an MPH pair route builds no matrix a compatible MPH analysis row has not already built. The MPH matrices are the same keyed artifacts the unary `mph_reml` and `mph_reml_ldms` routes use, so a relationship request and an analysis row that agree on cohort, genotype input and matrix settings share one built family; and at equal `ld_score_region_kb`, `ld_bins` and `ldms_maf_edges`, `mph_bivariate_reml_ldms` and `gcta_bivariate_reml_ldms` share one LD-by-MAF component plan and therefore one ordered SNP membership set. MPH matrices remain a different on-disk format from GCTA's and are never interchanged.

The sample is not the only thing that has to match before the two are the same model. MPH's REML is unconstrained — it returns a negative variance component where the likelihood puts one, and offers no option to constrain — while the pipeline's GCTA bivariate routes run GCTA's constrained default, because their rendered arguments are the endpoint-aware prevalence and your own `native_args` and nothing else. Wherever that constraint binds, and on a low-signal pair it does, the two routes are answering with different estimators rather than with two implementations of one. `--reml-no-constrain` is a legal named addition on a GCTA pair request and is the setting that makes the comparison matched; the pipeline does not add it for you, and adding it changes what that GCTA result is.

Treat both routes as a recommended REML candidate rather than an established default. Until parity with the GCTA pair routes has been demonstrated on data you recognise — which means a complete pair, where the two estimators fit the same individuals, and a matched constraint setting — read an MPH pair estimate as a second opinion rather than as the primary result.

A quantitative pair may request the REML and HE dense estimators together in one relationship row:

```csv title="relationship_manifest.csv"
relationship_id,left_analysis_id,right_analysis_id,left_summary_statistics_id,right_summary_statistics_id,relationship_methods,pair_quant_covariates,pair_cat_covariates
height_bmi,height,bmi,,,"gcta_bivariate_reml,gcta_bivariate_he",,
```

The pair phenotype is a deterministic full union of the two prepared endpoint sample sets in `FID`,`IID` order, with `NA` on a side where that trait is missing. It is prepared once per relationship and shared by every individual-level pair estimator that selects it, so selecting several of them defines one endpoint sample set and not one per method. What an estimator then does with a row missing one side is its own contract, and GCTA and MPH differ there, as described above. Pair covariates belong to the relationship, not either endpoint analysis. The primary pair request ID is deterministic:

```text
gcta_bivariate_reml--<relationship_id>
gcta_bivariate_reml_ldms--<relationship_id>
gcta_bivariate_he--<relationship_id>
gcta_bivariate_he_ldms--<relationship_id>
mph_bivariate_reml--<relationship_id>
mph_bivariate_reml_ldms--<relationship_id>
```

Summary pair request IDs use the same rule:

```text
<method>--<relationship_id>
```

### Method capabilities

Every selector carries only capabilities that validation, routing, reporting or GWASLab adaptation actively consumes. Selection stays explicit: the pipeline never substitutes one estimator for another according to sample size, memory, trait type or a failed task.

| Method group                  | Estimator family           | Input backend                        | Component model   | Trait support                         | Prevalence contract                   | Genotype bundle                       |
| ----------------------------- | -------------------------- | ------------------------------------ | ----------------- | ------------------------------------- | ------------------------------------- | ------------------------------------- |
| `regenie`                     | Whole-genome regression    | Direct PLINK genotypes (BED or PGEN) | None              | Quantitative and binary               | Not consumed                          | PLINK (BED or PGEN)                   |
| `gcta_fastgwa`                | Mixed linear model         | Sparse GRM                           | Single            | Quantitative and binary               | Not consumed                          | PLINK (BED or PGEN)                   |
| `ldak_kvik`                   | Mixed linear model         | Direct PLINK 1 genotypes             | Single            | Quantitative and binary               | Not consumed                          | PLINK 1                               |
| `gcta_greml`                  | REML                       | Dense GRM                            | Single            | Quantitative and binary               | Population value consumed             | PLINK (BED or PGEN)                   |
| `gcta_greml_ldms`             | REML                       | LDMS GRM family                      | LD/MAF stratified | Quantitative and binary               | Population value consumed             | PLINK 1                               |
| `ldak_reml`                   | REML                       | LDAK kinship                         | Single            | Quantitative and binary               | Population value consumed             | PLINK 1                               |
| `ldak_he`                     | Moment HE                  | LDAK kinship                         | Single            | Quantitative and binary               | Population value consumed             | PLINK 1                               |
| `ldak_pcgc`                   | PCGC                       | LDAK kinship                         | Single            | Binary only                           | Population value required             | PLINK 1                               |
| `ldak_fast_he`                | Moment HE (randomised)     | Direct PLINK 1 genotypes             | Single            | Quantitative and binary               | Population value consumed             | PLINK 1                               |
| `ldak_fast_pcgc`              | PCGC (randomised)          | Direct PLINK 1 genotypes             | Single            | Binary only                           | Population value required             | PLINK 1                               |
| `mph_reml`                    | REML (randomised trace)    | MPH GRM                              | Single            | Quantitative only                     | Not consumed                          | PLINK 1                               |
| `mph_reml_ldms`               | REML (randomised trace)    | MPH GRM family                       | LD/MAF stratified | Quantitative only                     | Not consumed                          | PLINK 1                               |
| GCTA bivariate REML           | REML                       | Dense or LDMS GRM                    | Single or LD/MAF  | Quantitative and binary               | Population value consumed             | PLINK (BED or PGEN); PLINK 1 for LDMS |
| GCTA bivariate HE             | Moment HE                  | Dense or LDMS GRM                    | Single or LD/MAF  | Quantitative only; no pair covariates | Not consumed                          | PLINK (BED or PGEN); PLINK 1 for LDMS |
| MPH bivariate REML            | REML (randomised trace)    | MPH GRM or GRM family                | Single or LD/MAF  | Quantitative only                     | Not consumed                          | PLINK 1                               |
| `ldsc_h2`, `ldsc_rg`          | LD-score regression        | Summary statistics                   | Single            | Quantitative and binary               | Population and sample values consumed | —                                     |
| `ldak_sumher`, `ldak_sumcors` | Summary tagging regression | Summary statistics                   | Tagging bundle    | Quantitative and binary               | Population and sample values consumed | —                                     |

The genotype bundle is the representation the estimator's own executable, or its matrix builder, reads, and it is what decides whether a cohort needs a PLINK 1 view derived for it: a run selecting only `PLINK (BED or PGEN)` methods never projects one, whatever representation the cohort was supplied in. `—` marks a summary estimator, which reads no genotypes at all.

Matrix-backed routes reuse a compatible matrix across requests and publish it only under `--save_relatedness_matrices`. A direct-genotype route requests no matrix at all: `ldak_fast_he` and `ldak_fast_pcgc` read the cohort's PLINK 1 view and build nothing under `quality_control/`.

A stochastic estimator randomises rather than enumerating, so two runs of the same input disagree unless the run is seeded. What varies differs by route. `ldak_fast_he` and `ldak_fast_pcgc` approximate the whole estimate, so both the heritability and its standard error move. `mph_reml`, `mph_reml_ldms` and the two MPH pair routes randomise only the trace terms of an otherwise exact REML fit, but that is enough to move their point estimate as well as its standard error. Each of these routes has a pipeline seed control — `ldak.fast_seed` and `mph.seed` — and leaving the LDAK one unset keeps the tool's own unseeded default. LDAK records no seed for an unseeded run, so the estimate cannot be reproduced afterwards. MPH's own unset default is the fixed seed `0`, which it echoes. For every MPH route a seed is necessary but not sufficient: see the option tables below. `ldak_he`, `ldak_pcgc` and `ldak_kvik` are stochastic only in their standard errors — their point estimates are reproducible — and the pipeline exposes no seed option for them yet, so repeated runs of those routes will report slightly different standard errors.

The component model names the variance-component structure the estimator fits, and is what decides whether a row or a pair request may configure the LD- and MAF-stratified plan settings.

`ldak_pcgc` and `ldak_fast_pcgc` are binary-only because LDAK requires a binary trait in those modes. `ldak_he` and `ldak_fast_he` accept quantitative and binary traits; a binary analysis with `population_prevalence` also receives native liability-scale estimates and conversion factors. Without prevalence, these HE routes report the observed scale. `mph_reml`, `mph_reml_ldms`, `mph_bivariate_reml` and `mph_bivariate_reml_ldms` are quantitative-only because MPH has no liability-scale conversion. For a binary or mixed pair, select `gcta_bivariate_reml` or `gcta_bivariate_reml_ldms`.

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

The [mixed summary-statistics manifest](../assets/examples/relational/summary_statistics_manifest.csv) illustrates both origins: one internal REGENIE result from `heterogeneous_qt` and external results loaded through the explicit `gwaslab` format. The companion [summary relationship](../assets/examples/relational/relationship_manifest_summary.csv), [request options](../assets/examples/relational/method_options_summary.json), and [reference-catalog shape](../assets/examples/relational/reference_catalog.json) show the complete declaration surface. Replace every `/refs/...` value in the catalog with a locally available scientific reference before launching; the pipeline deliberately rejects unavailable paths and does not infer a bundle from ancestry.

#### Common-variant meta-analysis

Declare a meta-analysis as another row in `--summary_statistics_manifest`. Its `summary_statistics_id` is
both the derived result identity and its request ID. Put at least two distinct declared base-summary IDs in
`source_summary_statistics_ids`, in the study order you want retained. Parents can be external summaries,
pipeline-generated summaries, or a mixture. Declare each parent in this manifest, even when an analysis
produces it. Meta-analysis-derived results cannot themselves be parents of another meta-analysis.

Leave `source`, `source_format`, `producer_analysis_id`, `producer_association_method`, `trait_id`,
`trait_type`, `genome_build`, `ancestry` and `source_method` blank on the derived row. Trait, type and build
must agree across its parents. The result inherits their common ancestry label, or `MULTI` when labels differ.
A multi-ancestry request explicitly selects `mrmega`; the labels describe the studies, while MR-MEGA derives
its continuous ancestry axes from the aligned per-study allele frequencies.

| Model    | Result                                                                           |
| -------- | -------------------------------------------------------------------------------- |
| `fixed`  | GWASLab inverse-variance pooled effect; supplies canonical `BETA`, `SE` and `P`. |
| `random` | Optional GWASLab conventional DerSimonian–Laird random-effects fields.           |
| `re2`    | Optional METASOFT Han–Eskin heterogeneous-effect association fields.             |
| `mrmega` | Optional MR-MEGA ancestry-aware meta-regression and heterogeneity fields.        |

A blank `meta_analysis_models` selects `fixed`. When you list models, include `fixed` alongside every optional
model. No model replaces another automatically. A blank `min_studies` means 2; otherwise supply an integer
from 2 through the declared parent count. The derived summary retains variants with at least this many
contributors to the native fixed-effect result. The threshold does not change native model files.

For example, add this row using the complete summary-manifest header:

```csv
summary_statistics_id,source_summary_statistics_ids,meta_analysis_models,min_studies
combined_eur,"cohort_a--regenie,published_eur","fixed,random,re2",2
```

The [complete two-study example](../assets/examples/relational/summary_statistics_manifest_meta.csv) uses small illustrative tables and can be launched from the repository root:

```bash
nextflow run . -profile docker \
    --summary_statistics_manifest assets/examples/relational/summary_statistics_manifest_meta.csv \
    --outdir results/meta_example
```

For an MR-MEGA request, declare the axis count explicitly in the fourth method-options namespace:

```json
{
  "meta_requests": {
    "combined_multi": {
      "mrmega": { "axes": 1 }
    }
  }
}
```

The pipeline forwards the integer `axes` value to MR-MEGA. Native model behavior remains native: METASOFT's
RE2 correction table covers 2–50 contributing studies and its very small p-values can be printed as zero;
MR-MEGA determines per-marker fit availability from the number of contributors and configured axes. Native
zeros, missing statistics and model-specific diagnostics are retained in the published results.

Same-ancestry derived summaries can select existing `heritability_methods` and be named as summary endpoints
in `relationship_manifest`. Current LDSC and LDAK summary methods require an ancestry-specific reference
strategy and therefore do not accept a `MULTI` derived result. Specify result-level binary prevalence when
needed by downstream methods; sample prevalence is not inferred or averaged from parents. Population
prevalence may inherit when all non-null parent declarations agree. Matching trait IDs assert compatible
phenotype meaning and effect scale; the pipeline does not infer sample independence from those labels.

Native results are published per selected model under `requests/<model>/<request_id>/`. The consolidated
fixed-effect summary is published under `summary_statistics/<summary_statistics_id>/` and contains the
selected optional model fields without replacing its fixed-effect `BETA`, `SE` or `P`.

## Advanced method options

`--method_options` is optional. The established form keeps its JSON root keyed by `analysis_id`; each value may contain `gcta`, `ldak`, `mph` and/or `regenie`. It remains supported unchanged. A namespaced document places those same entries under `analyses`, summary unary settings under `unary_requests`, relationship settings under `pair_requests`, and meta-analysis settings under `meta_requests`. Unlisted analysis settings receive their defaults. Every selected LDAK or LDSC summary request must explicitly choose a `reference_bundle_id`; nothing is inferred from the summary's ancestry label.

All four GCTA pair routes (`gcta_bivariate_reml`, `gcta_bivariate_reml_ldms`, `gcta_bivariate_he`, `gcta_bivariate_he_ldms`) expose `native_args` as an array of individual non-file GCTA tokens on the deterministic request ID. A relationship's deterministic LDMS request additionally owns its matrix construction settings; it never inherits them from either endpoint's unary analysis:

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

The wrapper rejects whitespace or shell syntax, path separators, environment assignments, undeclared file-like values, file-bearing invocation mechanics such as `--keep` and `--extract`, wrapper-owned flags such as `--grm`, `--pheno`, `--out`, `--reml-bivar`, `--reml-bivar-prevalence` and `--HEreg-bivar`, and flags selecting another primary GCTA operation such as `--pca`. Arguments remain native scientific options: the pipeline records them and presents all resulting estimates; it does not choose a preferred result.

`gcta_bivariate_reml_ldms` publishes GCTA's native per-stratum correlations and no genome-wide total. Your `--reml-no-constrain` and `--reml-maxit` choices remain yours, are rendered into the fit unchanged, and are recorded in `native.log`. On a small cohort a multi-stratum bivariate fit often reaches the constrain boundary and exhausts the iteration cap under GCTA's constrained default; `--reml-no-constrain` is what makes such a fit converge, and it changes what the result is.

`native_args` on a pair request is firewalled against the option namespace of the tool that actually runs, so the two families never share a list and a GCTA pair request in the same document keeps the GCTA firewall described above. An MPH pair request is admitted by a positive allow-list instead: only `--min_maf`, `--min_hwe_pval` and `--verbose` may be set. Every other MPH option — the invocation mode, the matrix list, the trait and covariate name lists, the output basename, the thread count, and every curated stochastic and solver control — is rendered by the wrapper from validated request state and is refused at ingress, as is an abbreviation of an accepted option and any option MPH does not have.

`mph_bivariate_reml` and `mph_bivariate_reml_ldms` additionally accept the curated `iterations`, `tolerance`, `random_vectors`, `seed` and `save_memory` settings directly on the pair request, under the same names, types and rules as the analysis-namespace `mph` family documented below. `mph_bivariate_reml_ldms` also owns its `ld_score_region_kb`, `ld_bins` and `ldms_maf_edges` matrix settings, and inherits them from neither endpoint's unary analysis:

```json
{
  "pair_requests": {
    "mph_bivariate_reml--height_bmi": {
      "random_vectors": 1000,
      "seed": 7
    },
    "mph_bivariate_reml_ldms--height_bmi": {
      "ld_score_region_kb": 200,
      "ld_bins": 4,
      "ldms_maf_edges": [0, 0.01, 0.05, 0.2, 0.5],
      "random_vectors": 1000,
      "seed": 7,
      "native_args": ["--min_maf", "0.01"]
    }
  }
}
```

Summary requests use the same deterministic ownership boundary. LDAK receives exactly one staged `tagging_file`; LDSC receives separate staged `hapmap3_snplist`, `reference_ld_scores` and `regression_weights` roles, the latter two declared as LDSC stems rather than directories (see [Reference catalog](#reference-catalog)). `native_args` may contain non-file scientific tokens only. Wrapper-owned operation, input, output and thread flags are rejected, as are every LDSC option that selects an alternate operation or consumes an undeclared file role. These structural rejections apply to both bare `--option value` and inline `--option=value` forms without depending on whether a path exists or resembles a known extension.

For `ldak_sumher` and `ldak_sumcors`, the pipeline adapts each distinct GWASLab summary once to LDAK's `Predictor A1 A2 Z n A1Freq` contract, with `A1` equal to the GWASLab effect allele and `Z = BETA / SE`; the GWASLab artifact remains unchanged. Both routes use `--cutoff 0.01` unless a request explicitly supplies `--cutoff` or `--truncate`. SumCors initially accepts `LDAK-Thin`, `Uniform-GCTA` and `Human-Default` tagging bundles. Binary SumHer receives population prevalence and sample ascertainment only when both are declared. SumCors receives the two ordered prevalence/ascertainment pairs only when both endpoints are binary and all four values are present; mixed-trait and incomplete binary pairs run without liability arguments. LDAK's native ambiguous-variant exclusion and complete-summary checks remain enabled unless an accepted scientific override changes them.

Each LDSC H2 request runs once and publishes one native log. A binary summary receives the prevalence flags for native liability conversion only when both `sample_prevalence` and `population_prevalence` are declared; otherwise it runs on the observed scale. LDSC RG also runs once per ordered pair, supplying prevalence only when at least one endpoint is binary and every binary endpoint declares both values. Quantitative endpoints in that mixed invocation use LDSC's native `nan` prevalence placeholder. The pipeline does not parse these logs into a common heritability, covariance or correlation family. Request-level `native_args` cannot override either prevalence flag.

```json
{
  "unary_requests": {
    "ldsc_h2--height--regenie": {
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
    "ldsc_h2--height--regenie": {
      "reference_bundle_id": "ldsc_eur"
    },
    "ldsc_h2--height--regenie--alternate": {
      "primary_request_id": "ldsc_h2--height--regenie",
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
    },
    "ldsc_eur_baselineld": {
      "genome_build": "GRCh37",
      "ancestry": "EUR",
      "variant_id_system": "rsid",
      "hapmap3_snplist": "/refs/w_hm3.snplist",
      "reference_ld_scores": "/refs/baseline/baselineLD.",
      "regression_weights": "/refs/weights/weights.hm3_noMHC."
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

`reference_ld_scores` and `regression_weights` are LDSC stems, written exactly as they would be typed after `--ref-ld-chr` and `--w-ld-chr`. LDSC completes a stem with `<chr>.l2.ldscore.gz`, and the reference set additionally with `<chr>.l2.M_5_50`. The pipeline splits the stem at the last `/`: the text before it is the directory it stages, and the text after it is the file-name prefix it hands back to LDSC. `/refs/eur_w_ld_chr/` therefore names `/refs/eur_w_ld_chr/1.l2.ldscore.gz` with an empty prefix, exactly as it always has, and `/refs/baseline/baselineLD.` names `/refs/baseline/baselineLD.1.l2.ldscore.gz`. A bare directory path with no trailing `/` is read as a stem too, so `/refs/eur_w_ld_chr` would look for `/refs/eur_w_ld_chr1.l2.ldscore.gz`; write the trailing `/` when the files are named `<chr>.l2.*`. Nothing is inferred from the directory listing: the declared text is authoritative.

The catalog may also declare a SHA-256 digest beside each role. Preflight checks the document shape, family, required roles, digest syntax and path availability, and for each LDSC stem that its directory exists and holds the chromosome-1 files the stem resolves to. In the first release, `genome_build`, `ancestry`, `variant_id_system` and `model` are recorded request metadata rather than a pipeline certification of scientific compatibility. The user owns reference selection.

| GCTA option          | Type and default                      | Consumer and constraints                                                                                |
| -------------------- | ------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| `grm_maf`            | Number or `null`; `null`              | Dense `gcta_greml` matrix filter; `0` to `0.5` inclusive when set.                                      |
| `grm_extract`        | Resource path or absent; absent       | Predictor list staged for dense `gcta_greml` matrix construction.                                       |
| `reml_no_constrain`  | Boolean; `false`                      | `gcta_greml` and `gcta_greml_ldms`; disables variance-component constraints.                            |
| `sparse_cutoff`      | Number; `0.05`                        | `gcta_fastgwa` only; `0` to `1` inclusive.                                                              |
| `ld_score_region_kb` | Positive integer; `200`               | `gcta_greml_ldms` only; LD-score window in kilobases.                                                   |
| `ld_bins`            | Positive integer; `4`                 | `gcta_greml_ldms` only; number of LD-score strata.                                                      |
| `ldms_maf_edges`     | Number array; `[0,0.01,0.05,0.2,0.5]` | `gcta_greml_ldms` only; strictly increasing, at least two values, beginning at `0` and ending at `0.5`. |

The same three LDMS setting names are accepted under a `gcta_bivariate_reml_ldms--<relationship_id>`, `gcta_bivariate_he_ldms--<relationship_id>` or `mph_bivariate_reml_ldms--<relationship_id>` pair request. Their defaults are `200`, `4` and `[0,0.01,0.05,0.2,0.5]`. The resolved values become part of the matrix reuse key, so a unary and pair request reuse one matrix family only when cohort, genotype input and all three settings agree exactly — independent of which LDMS pair method token requested it.

Those three settings define one LD-by-MAF component plan: the LD scores of a cohort's genotype view and the ordered, disjoint SNP groups derived from them. The plan is its own reusable artifact, keyed by the genotype view and the three settings alone, so the GCTA and the MPH stratified families that declare one set of settings share one plan and therefore one SNP membership set instead of partitioning the same variants twice, whether they were requested by an analysis row, by a relationship, or by both. On an analysis row the settings stay in the `gcta` family because GCTA's own LD-score pass builds the plan: a row selecting only `mph_reml_ldms` still declares `ld_score_region_kb`, `ld_bins` and `ldms_maf_edges` under `gcta`. A relationship instead declares them on its own LDMS pair request, as described above. With `--save_relatedness_matrices` each plan is published once under `quality_control/ldms_component_plans/<plan_key>/`, separately from the matrix families that consume it.

| LDAK option          | Type and default                | Consumer and constraints                                                                                                                                                                                                                                                                                                                                                |
| -------------------- | ------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `model`              | String; `human_default`         | Kinship and direct-genotype routes; `human_default` or `custom`.                                                                                                                                                                                                                                                                                                        |
| `power`              | Number; `-0.25`                 | Kinship and direct-genotype routes; `-2` to `0`. `human_default` fixes `-0.25`; select `custom` for another value.                                                                                                                                                                                                                                                      |
| `weights_policy`     | String; `equal`                 | Kinship and direct-genotype heritability routes; `equal` applies no predictor weights, while `provided` requires the staged `weights` resource.                                                                                                                                                                                                                         |
| `weights`            | Resource path or absent; absent | Kinship and direct-genotype routes; accepted exactly with `weights_policy: provided`. LDAK gives weight zero to any predictor absent from the file and only warns in the native `.log` (`contains weights for only`), so supply a weight for every predictor of the cohort.                                                                                             |
| `relatedness_filter` | Boolean; `false`                | Kinship routes only; optionally derives an unrelated subset. Rejected on a row that selects no kinship method. On a mixed row the pipeline warns: the direct-genotype estimators cannot consume the keep list, so the two estimates use different sample sets.                                                                                                          |
| `fast_repetitions`   | Positive integer or `null`      | `ldak_fast_he` and `ldak_fast_pcgc` only; `null` keeps LDAK's own 100 random vectors.                                                                                                                                                                                                                                                                                   |
| `fast_num_blocks`    | Integer or `null`               | `ldak_fast_he` and `ldak_fast_pcgc` only; `null` keeps LDAK's own 200 blocks. These are **predictor** jackknife blocks — LDAK partitions predictors, not samples — so the setting moves the standard error and not the point estimate. A block-jackknife standard error is not comparable across block settings, nor to the exact `ldak_he`/`ldak_pcgc` standard error. |
| `fast_seed`          | Integer or `null`               | `ldak_fast_he` and `ldak_fast_pcgc` only; `null` leaves the run unseeded. An unseeded LDAK run records no seed anywhere, so its estimate cannot be reproduced afterwards.                                                                                                                                                                                               |
| `kvik_step1_subset`  | String; `all`                   | `ldak_kvik` only; `all`, `thin_common`, or `provided`.                                                                                                                                                                                                                                                                                                                  |
| `predictor_extract`  | Resource path or absent; absent | `ldak_kvik` only; required exactly with `kvik_step1_subset: provided`.                                                                                                                                                                                                                                                                                                  |
| `kvik_step2_keep`    | Resource path or absent; absent | `ldak_kvik` Step 2 only; optional FID/IID sample keep file passed to native `--keep`.                                                                                                                                                                                                                                                                                   |

Every LDAK analysis-row estimator — `ldak_kvik`, `ldak_reml`, `ldak_he`, `ldak_pcgc`, `ldak_fast_he` and `ldak_fast_pcgc` — requires a complete covariate cell for every phenotyped sample. LDAK never parses a missing `--covar` cell: it leaves whatever value was last in its read buffer in place, so the fitted covariate silently becomes the neighbouring column's value for that sample, the previous row's last value, or uninitialised memory. A missing `--factors` cell instead becomes an additional factor level. In neither case does LDAK drop the sample or say so. The pipeline therefore rejects such an analysis at phenotype preparation and names the offending `FID IID column` cells, rather than reporting a covariate model nobody declared.

| MPH option       | Type and default                  | Consumer and constraints                                                                                                                                                                                                                                                                                                       |
| ---------------- | --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `iterations`     | Positive integer or `null`        | `mph_reml` and `mph_reml_ldms` only; `null` keeps MPH's own 20 REML iterations. Exhausting them is a `Warning: not converged after N iterations.` line at exit 0 with a complete, well-formed and materially different result written beside it; the run is not failed on it, so check the native log before using the result. |
| `tolerance`      | Number greater than `0` or `null` | `mph_reml` and `mph_reml_ldms` only; `null` keeps MPH's own `0.01` stopping rule on the predicted log-likelihood change. The native log records the options and the solver trace records completed iterations.                                                                                                                 |
| `random_vectors` | Positive integer or `null`        | `mph_reml` and `mph_reml_ldms` only; `null` keeps MPH's own 100 random vectors. This is the one setting that narrows the seed-to-seed spread of the estimate; see the paragraph below.                                                                                                                                         |
| `seed`           | Integer or `null`                 | `mph_reml` and `mph_reml_ldms` only; `null` keeps MPH's own default, which is the fixed integer `0` rather than entropy. The native log records the effective value.                                                                                                                                                           |
| `save_memory`    | Boolean; `false`                  | `mph_reml` and `mph_reml_ldms` only; passes `--save_memory`, which reads the relationship matrices from disk during the fit instead of holding them in memory. Measured with covariates it left `var`, `seV` and `pve` identical to every printed digit and moved `enrichment` at the sixth significant figure.                |

That table is the analysis-row `mph` family. The same five settings are accepted with identical names, types and rules on an `mph_bivariate_reml--<relationship_id>` or `mph_bivariate_reml_ldms--<relationship_id>` pair request, where they are declared directly on the request rather than under a family key.

An unseeded MPH run is reproducible, which is not true of the other stochastic heritability routes: MPH's unset seed is the fixed integer `0`, not entropy, so repeating a run on one machine at one thread count returns the same estimate. What a seed does not buy is stability. Any value moves the point estimate itself through the stochastic trace, and the reported standard error carries none of that movement: measured on a 200-sample cohort over seeds 1 to 8, the proportion of variance explained spanned `0.114` at 50 random vectors (`0.022` to `0.137`), `0.052` at MPH's default of 100 (`0.064` to `0.116`) and `0.016` at 1000 (`0.064` to `0.079`), while the reported standard error moved only between `0.124` and `0.129` at that default. Raise `mph.random_vectors` for anything you intend to publish. Two further things move the estimate, both in the sixth to seventh significant digit and neither of them a method option: the thread count the executor gives the task, and `mph.save_memory`. The native log records the effective settings, so a published MPH estimate is reproducible against the same seed, random-vector count, thread count and memory mode, and not otherwise.

MPH builds its own relationship matrices and takes no construction option for them: no MAF filter, no predictor extract list and no scaling exponent enters one, and MPH's only native quality filter is to drop a variant whose minor allele frequency is zero. The variant universe is instead declared: a variant outside the autosomes enters the matrix with weight zero, which MPH treats exactly as an omitted one. `gcta.grm_maf` and `gcta.grm_extract` restrict the GCTA dense matrix and nothing else, so a row that supplies either alongside an MPH selector is refused rather than routed; otherwise the run would publish a filtered GCTA estimate and an unfiltered MPH one side by side under one analysis identifier, with nothing in either native output saying so.

MPH's covariate interface shapes what the pipeline hands it in two ways. MPH synthesises an intercept only when no covariate is named at all, so the pipeline always writes an explicit all-ones `intercept` column and names it first; and MPH never expands a categorical covariate, so factor columns are dummy-encoded by the same rule the LDAK matrix-adjustment design uses. When the resulting design is rank deficient MPH prints one warning, prunes columns until it is not, and names none of the columns it dropped; measured, one collinear design dropped that explicit `intercept`, turning a covariate-adjusted fit into a fit through the origin. That is a different model rather than a rounding difference: fitting the same two covariates with and without an intercept moved the proportion of variance explained from `0.0849` to `0.1083`. Inspect MPH's native fixed-effect result and log for the covariates retained in the fit.

Unlike every LDAK analysis-row estimator, MPH does not need a complete covariate cell. An empty field is its only missing representation, and a sample carrying one is dropped: MPH reports the reduced set as `Non-missing analysis set contains N individuals` in its native log. Preparation rewrites the pipeline's missing spellings to empty fields.

Both MPH routes run on a `linux/amd64` image and only there: the binary is statically linked against Intel oneMKL, there is no Conda package and no build for another architecture, so a task scheduled elsewhere fails at exec rather than falling back to one.

| REGENIE option      | Type and default         | Consumer and constraints                                                                             |
| ------------------- | ------------------------ | ---------------------------------------------------------------------------------------------------- |
| `step1_bsize`       | Positive integer; `1000` | Block size for the analysis's own Step 1 fit.                                                        |
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

Resource paths are resolved and staged at ingress. Step 1 resources contribute to the matrix, predictor, or fit identity that consumes them; `ldak.kvik_step2_keep` is a Step 2-only input and cannot fragment preparation or Step 1 reuse. `gcta_grm_parts` is operational partitioning and remains configuration/profile-only, never a method option.

The authoritative structural contracts are [`schema_cohort_manifest.json`](../assets/schema_cohort_manifest.json), [`schema_analysis_manifest.json`](../assets/schema_analysis_manifest.json), [`schema_summary_statistics_manifest.json`](../assets/schema_summary_statistics_manifest.json), [`schema_relationship_manifest.json`](../assets/schema_relationship_manifest.json), and [`schema_reference_catalog.json`](../assets/schema_reference_catalog.json); origin-, relationship-, request-, method-, trait- and resource-aware diagnostics come from the central preflight validator.

### Validation diagnostics

Structural failures name the manifest and invalid column. Cross-row preflight failures additionally name the CSV row and `cohort_id`, `analysis_id`, `summary_statistics_id` or `relationship_id`: examples include an incomplete genotype group, conflicting duplicates, an invalid or mixed summary origin, an internal producer mismatch, orphan endpoints, same-endpoint or same-trait pairs, cross-cohort GCTA pairs, reversed duplicates, unknown or repeated method tokens, invalid binary coding, route-inapplicable prevalence, a binary or mixed pair selecting a quantitative-only method such as `gcta_bivariate_he`, `gcta_bivariate_he_ldms`, `mph_bivariate_reml` or `mph_bivariate_reml_ldms`, and declared pair covariates against a method with no native covariate parameter. Method-options and reference-catalog failures name the document, request or bundle ID, qualified option or resource role, and reason before task submission.

### Phenotype input preparation

`PREPARE_PHENOTYPE_INPUTS` prepares each selected phenotype in a shared two-identifier-plus-trait layout. Binary source values matching `control_value` and `case_value` become `0` and `1`; missing values and unmatched binary values become `NA`. Quantitative values must be numeric. Quantitative and categorical covariates remain distinct for tools with separate native interfaces.

Categorical covariates are passed to REGENIE by their names in the `cat_covariates` header and dummy-coded by REGENIE.

### Run-level defaults

Parameters divide by what they control. Scientific specification — what is estimated and under which
model — is declared per analysis in the method-options document, not on the command line. Execution
strategy — task partitioning, execution modes, job counts, block sizes and low-memory modes — is set by
the parameters in the schema's `Execution options` group, whose defaults live in `conf/execution.config`. CPU, memory and
time requests are not parameters at all: change them with a custom process configuration passed via `-c`.

| Parameter or behaviour            | Default and rationale                                                                                                                                                                                            |
| --------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--regenie_step1_mode`            | `standard`, the simplest one-task Step 1. Use `chunked` with `--regenie_step1_jobs` when a large cohort needs REGENIE's split-L0/run-L0/run-L1 execution family.                                                 |
| `--regenie_lowmem`                | `true`, keeping Step 1's temporary prediction blocks in the task work directory to reduce memory use.                                                                                                            |
| REGENIE scientific method options | Per-analysis `regenie.*` defaults enable approximate Firth fallback below `0.01` for binary traits and leave `min_mac` unset so REGENIE's own versioned policy applies.                                          |
| GWASLab reference parameters      | Unset. Every association output is still standardised; reference-dependent allele checks, rsID assignment and strand inference run only when you provide the corresponding build-specific FASTA or VCF resource. |
| Save controls                     | Off. Intermediates stay out of the results directory unless explicitly requested, avoiding unexpectedly large published output.                                                                                  |

The two opt-in save controls are:

- `--save_normalised_phenotypes`: publish headered prepared phenotype and covariate files under `phenotypes/<analysis_id>/`. The parameter name is retained for compatibility.
- `--save_relatedness_matrices`: publish each reusable GCTA, LDAK or MPH base or derived matrix artifact once under `quality_control/relatedness_matrices/<key>/`, and each LD-by-MAF component plan once under `quality_control/ldms_component_plans/<plan_key>/`.

REGENIE and LDAK-KVIK Step 1 bundles remain internal work outputs consumed directly by Step 2. Selective `-resume` requires preserving the Nextflow cache and work directory; the pipeline does not publish or import Step 1 bundles.

See the [output documentation](output.md) for the exact files and publication exceptions.

REGENIE runs Step 2 separately for each analysis. Analyses with identical genotype views, trait type,
phenotype and covariate bytes, and Step 1 block size share one Step 1 prediction bundle. A Step 2 option
change leaves that fit reusable; a phenotype change creates a new fit for the affected analysis.

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

| What you see                                              | Meaning                                                                                           | First evidence to inspect                                                                                                |
| --------------------------------------------------------- | ------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| An error naming a manifest, CSV row, identity and field   | Preflight rejected input before task submission                                                   | Read the complete grouped error report, fix every listed row and relaunch.                                               |
| `Method-options document ... request_id ... option ...`   | A deterministic or named request is incomplete, invalid or conflicts with the invocation firewall | Read the request namespace, request ID, option and reason; then correct the complete request configuration.              |
| `Reference catalog ... reference_bundle_id ... field ...` | A bundle is malformed, duplicated across families or has an unavailable required resource         | Correct the named bundle role; compatibility metadata is recorded but not scientifically certified by the pipeline.      |
| `Process ... terminated with an error`                    | A task was submitted and failed                                                                   | Inspect `.nextflow.log` and the task's `.command.err`, `.command.out` and `.command.log` in the reported work directory. |

### Manifest validation failed

Read the reported CSV path, row number, entity ID, field name and reason from left to right. The validator reports all linked-manifest errors it can find in one launch, so correct every bullet before rerunning.

| Diagnostic fragment                                                                                                       | Supported action                                                                                                                                                                                                                                                                         |
| ------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Validation of samplesheet failed` and a named field                                                                      | Correct the schema violation. Required columns must be present; schema-optional columns may be omitted; unexpected columns are rejected. Column order is not significant.                                                                                                                |
| `header row 1 repeats column ...`                                                                                         | Keep each column name once. Quoted header names are accepted.                                                                                                                                                                                                                            |
| `no genotype group is populated`, `a second genotype group is populated` or `genotype group ... is only partly populated` | Populate exactly one complete representation: `pgen`/`psam`/`pvar` or `bed`/`bim`/`fam`.                                                                                                                                                                                                 |
| `duplicate cohort_id` or `duplicate analysis_id`                                                                          | Keep one row for each identity. Conflicting duplicates also report the fields that differ.                                                                                                                                                                                               |
| `undefined cohort_id`                                                                                                     | Make the analysis row's `cohort_id` match one cohort-manifest identity exactly.                                                                                                                                                                                                          |
| `unknown method`, `listed more than once` or `row selects no method`                                                      | Use each documented selector at most once and populate at least one of `association_methods` or `heritability_methods`. See [Analysis manifest fields](#analysis-manifest-fields).                                                                                                       |
| A binary `case_value` or `control_value` diagnostic                                                                       | Supply both distinct source codes for a binary trait. Remove both from a quantitative row. The pipeline does not infer `trait_type` from phenotype values.                                                                                                                               |
| A `population_prevalence` diagnostic                                                                                      | Use it only for a binary heritability analysis whose selected estimator consumes it, and provide it for `ldak_pcgc` and `ldak_fast_pcgc`. See [Analysis manifest fields](#analysis-manifest-fields).                                                                                     |
| `support quantitative traits only` or `support binary traits only`                                                        | Match the row's `trait_type` to the selected estimator. The diagnostic lists the estimators that do support the declared trait type.                                                                                                                                                     |
| `missing covariate cells`                                                                                                 | Complete the covariate cells the diagnostic names, or drop those samples from the phenotype file. Every LDAK analysis-row estimator reads a missing cell as a value rather than excluding the sample.                                                                                    |
| `restricts the GCTA dense predictor set only`                                                                             | Remove `gcta.grm_maf` or `gcta.grm_extract`, or drop the MPH selector the diagnostic names from that row. Either option filters the GCTA dense matrix alone; MPH's matrices are built over the declared MPH variant universe, so the two estimates would cover different predictor sets. |
| `names repeat`                                                                                                            | Give every declared weight column a unique name.                                                                                                                                                                                                                                         |
| A summary origin or deterministic identity diagnostic                                                                     | Populate exactly one complete external, producer or meta-analysis origin. Internal IDs must be `<producer_analysis_id>--<producer_association_method>`. See [Summary-statistics manifest fields](#summary-statistics-manifest-fields).                                                   |
| An undefined or self-paired summary endpoint                                                                              | Declare each referenced summary ID, use distinct endpoint and trait IDs on the two sides, and ensure any combined analysis/summary side has recorded producer correspondence.                                                                                                            |

### Method-options validation failed

The diagnostic names the document, entity or request ID, fully qualified option and reason. Fix malformed JSON; use the `gcta`, `ldak`, `mph` and `regenie` analysis families; and use only the `analyses`, `unary_requests`, `pair_requests` and `meta_requests` namespaces described under [Advanced method options](#advanced-method-options). A request may configure only a method already selected by its summary or relationship declaration. Named additions must provide their own complete settings, including a reference bundle for LDAK or LDSC.

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
- Prepared genotypes, prepared phenotypes and covariates, and relatedness matrices are unpublished by default; enable the corresponding save control before expecting those directories. REGENIE and LDAK-KVIK Step 1 predictions are internal work outputs and have no publication control.
- GWASLab reference parameters are optional. Standardised association output is still produced without them, but reference-dependent allele checks, flips, rsID assignment and strand inference are not. `genome_build`, not `ancestry`, selects the build-specific resources.
- Published intermediates are retention outputs, not importable cross-run caches. Supported reuse requires retained Nextflow work and `-resume`.

### `-resume` did not reuse work

Nextflow reuses a task only when its inputs, pipeline code and relevant configuration still match a retained cache entry. Keep the original `work/` directory and `.nextflow/` cache, rerun from the same launch context where practical, and use `nextflow log` to find a run name for `-resume <run-name>`. A changed input file, manifest, method option, pipeline revision, process configuration or missing work directory can require recomputation. See [`-resume`](#-resume).

### Get help

Before requesting help, retain the pipeline version, exact launch command with secrets removed, `.nextflow.log`, the failing process name and relevant `.command.*` files, the reported manifest row or a minimal redacted reproducer, executor/profile and container runtime, and whether the original work directory remains available for `-resume`.

Use the [nf-core troubleshooting guide](https://nf-co.re/docs/running/troubleshooting), open a [GitHub issue](https://github.com/nf-core/gwas/issues) with a reproducible report, or ask in the nf-core [`#gwas` Slack channel](https://nfcore.slack.com/channels/gwas).
