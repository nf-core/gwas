# nf-core/gwas: Output

## Introduction

This document describes the files that nf-core/gwas publishes beneath `--outdir`. Native association, heritability and declared pairwise results are retained. Every internal association result and every external summary source converges on one canonical summary-statistics contract with a provenance sidecar. Run-level provenance is collected in MultiQC and `pipeline_info/`.

Intermediates are unpublished by default. The optional directories described below appear only when their corresponding save control is enabled.

### Naming and provenance

The shared result prefix grammar is:

```text
<analysis_id>.<method>[.<shard>]
```

`<analysis_id>` is copied from the analysis manifest and identifies one cohort-trait analysis unit. `<method>` is one of the method-selector tokens documented in [Usage](usage.md#relational-manifest-input). A producing tool can add a native result suffix after that prefix.

Summary results use a separate first-class identity. Pipeline-generated association summaries use `<analysis_id>--<association_method>`; external summaries use the declared `summary_statistics_id`. Both publish as `<summary_statistics_id>.canonical.tsv.gz` with `<summary_statistics_id>.provenance.json` in the same identity-addressed directory.

Use the following provenance chain for any result:

1. For an analysis result, read `<analysis_id>` and `<method>` from its parent directories and filename. Find that analysis row, follow its `cohort_id`, and inspect its method options.
2. For a summary result, read `summary_statistics_id` from its directory and provenance sidecar. The sidecar identifies an internal producer analysis/method or the external source basename and checksum, plus the transformation path.
3. Map the method to its producing tool using the table below.
4. Read tool versions from `pipeline_info/nf_core_gwas_software_mqc_versions.yml`. The pipeline version and complete run parameters are recorded by the `pipeline_info/` reports and `params_<timestamp>.json`.

Pairwise outputs instead use the deterministic request ID `<method>--<relationship_id>`. Find `relationship_id` in `--relationship_manifest`, follow its ordered left and right analysis IDs into `--analysis_manifest`, and use `requests/<method>/<request_id>/provenance.json` for the exact endpoint orientation, dense or LDMS matrix reuse key and native basename, effective prevalence, native arguments, all parsed native components, warnings and completion classification.

| Method token                                                                                       | Producing tool |
| -------------------------------------------------------------------------------------------------- | -------------- |
| `plink2`                                                                                           | PLINK 2        |
| `regenie`                                                                                          | REGENIE        |
| `gcta_fastgwa`, `gcta_greml`, `gcta_greml_ldms`, `gcta_bivariate_reml`, `gcta_bivariate_reml_ldms` | GCTA           |
| `ldak_kvik`, `ldak_reml`, `ldak_he`, `ldak_pcgc`                                                   | LDAK 6         |
| `ldak_sumher`, `ldak_sumcors`                                                                      | LDAK 6.3       |
| `ldsc_h2`, `ldsc_rg`                                                                               | LDSC           |

Together, the result prefix, retained cohort and analysis manifests, optional method-options document, and `pipeline_info/` artifacts identify the analysis, cohort, trait, genome build, method, scientific settings, pipeline version and producing tool version. Preserve them with an archived result.

Analysis attribution applies under `association/` and `heritability/individual/`; summary attribution applies under `summary_statistics/` and summary-level request outputs. Optional prepared genotypes and relatedness matrices are deliberately shared artifacts rather than trait-method results: `genotypes/` is attributed to `cohort_id`, while `quality_control/relatedness_matrices/` is attributed to its reuse key and may serve several analysis rows.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and publishes:

- [Association](#association)
  - [PLINK 2](#plink-2)
  - [REGENIE](#regenie)
  - [GCTA fastGWA-MLM](#gcta-fastgwa-mlm)
  - [LDAK-KVIK](#ldak-kvik)
- [Summary statistics](#summary-statistics)
- [Heritability](#heritability)
  - [GCTA GREML and GREML-LDMS](#gcta-greml-and-greml-ldms)
  - [LDAK estimators](#ldak-estimators)
- [Pairwise GCTA bivariate REML](#pairwise-gcta-bivariate-reml)
- [Quality control and optional prepared data](#quality-control-and-optional-prepared-data)
- [MultiQC](#multiqc)
- [Pipeline information](#pipeline-information)

## Association

Native association output is published by method and then analysis. The native tables are not rewritten; deterministic filename tokens added by tools are removed where necessary so every published name follows the shared prefix grammar.

### PLINK 2

<details markdown="1">
<summary>Output files</summary>

[PLINK 2](https://www.cog-genomics.org/plink/2.0/assoc) receives the normalised phenotype, preserves its native result columns and is configured to include `A1_FREQ`, `OBS_CT`, `BETA`, `SE` and `P` for harmonisation. The pipeline removes PLINK 2's fixed `.PHENO` token from the published filename; file content is unchanged.

- `association/plink2/<analysis_id>/`
  - `<analysis_id>.plink2.glm.linear`: Native PLINK 2 `--glm` result for a quantitative trait.
  - `<analysis_id>.plink2.glm.logistic.hybrid`: Native PLINK 2 `--glm` result for a binary trait with Firth fallback available.

</details>

For binary traits, Firth fallback is enabled by default so complete or quasi-complete separation, often encountered for rare variants, can still produce an estimate; this is why the result uses the `.logistic.hybrid` extension. The common normalised binary coding is `0`/`1`/`NA`, so the pipeline passes `--1`. It also uses `--covar-variance-standardize` when covariates are present to avoid PLINK 2's numerical-stability failure for differently scaled covariates; this invertible covariate reparameterisation does not change the reported genotype effect.

### REGENIE

<details markdown="1">
<summary>Output files</summary>

[REGENIE](https://rgcgithub.github.io/regenie/) fits a whole-genome prediction model in Step 1 and tests variants in Step 2. Standard Step 1 is the default because it is the simplest execution path; `--regenie_step1_mode chunked` and `--regenie_step1_jobs` use split-L0/run-L0/run-L1 when a large cohort needs work divided into smaller jobs. `--regenie_lowmem` defaults to `true` so temporary prediction blocks stay in the task work directory rather than memory.

- `association/regenie/<analysis_id>/`
  - `<analysis_id>.regenie.gz`: Native, space-delimited REGENIE Step 2 association result.
- `intermediates/association_predictions/regenie/<analysis_id>/` (with `--save_association_predictions`)
  - `<analysis_id>.regenie_step1_pred.list`: Prediction-list manifest from standard Step 1 or chunked Run L1.
  - `<analysis_id>.regenie_step1_<chromosome>.loco.gz`: Leave-one-chromosome-out prediction table.

</details>

For binary traits, the per-analysis `regenie.firth`, `regenie.firth_approx` and `regenie.firth_p_threshold` method options default to approximate Firth correction below `0.01`. `regenie.min_mac` defaults to `null`, retaining REGENIE's own versioned minimum-MAC policy unless an analysis explicitly overrides it.

The published Step 2 file is native output with the fixed `_PHENO` token removed from its name. Step 1 predictions are reusable intermediates and remain unpublished unless `--save_association_predictions` is enabled; chunk-planning files, temporary low-memory predictions and logs remain in the work directory.

### GCTA fastGWA-MLM

<details markdown="1">
<summary>Output files</summary>

[GCTA](https://yanglab.westlake.edu.cn/software/gcta/) runs `--fastGWA-mlm` for quantitative traits and `--fastGWA-mlm-binary` for binary traits. Plain fastGWA linear regression is not selectable. The method option `gcta.sparse_cutoff` defaults to `0.05`, following the official fastGWA example, and controls construction of its sparse relatedness matrix. The native table is passed unchanged to GWASLab.

- `association/gcta_fastgwa/<analysis_id>/`
  - `<analysis_id>.gcta_fastgwa.fastGWA`: Native GCTA fastGWA-MLM result.

</details>

### LDAK-KVIK

<details markdown="1">
<summary>Output files</summary>

[LDAK-KVIK](https://dougspeed.com/ldak-kvik/) fits a Step 1 prediction model from the PLINK 1 compatibility bundle prepared once per cohort and tests the full bundle in Step 2. `ldak.kvik_step1_subset` defaults to `all`; `thin_common` requests deterministic thinning and `provided` requires the stageable `ldak.predictor_extract` resource. The choice changes prediction reuse identity. The native `.assoc` file is published unchanged. The reusable three-file Step 1 bundle is optional output; effects, thinning progress and logs remain in the work directory.

- `association/ldak_kvik/<analysis_id>/`
  - `<analysis_id>.ldak_kvik.step2.assoc`: Native LDAK-KVIK Step 2 association table.
- `intermediates/association_predictions/ldak_kvik/<analysis_id>/` (with `--save_association_predictions`)
  - `<analysis_id>.ldak_kvik.step1.root`: LDAK Step 1 model root file, attributed to the requesting analysis.
  - `<analysis_id>.ldak_kvik.step1.loco.details`: Leave-one-chromosome-out model details used by LDAK-KVIK Step 2.
  - `<analysis_id>.ldak_kvik.step1.loco.prs`: Leave-one-chromosome-out polygenic scores used by LDAK-KVIK Step 2.

</details>

## Summary statistics

<details markdown="1">
<summary>Output files</summary>

[GWASLab](https://cloufield.github.io/gwaslab/) standardises every pipeline-generated association result and every external source declared with `source_mode: raw`; native association results remain available under `association/`. An external `source_mode: canonical` table bypasses GWASLab but passes through the same canonical validator. The external source itself and the temporary GWASLab table are not republished, so each scientific summary result appears only once.

- `summary_statistics/<summary_statistics_id>/`
  - `<summary_statistics_id>.canonical.tsv.gz`: Gzip-compressed, tab-delimited `nfcore_gwas_canonical_v1` table.
  - `<summary_statistics_id>.provenance.json`: Safe source, producer, transformation, checksum and canonical-contract provenance.

</details>

The required columns are `SNPID`, `CHR`, `POS`, `EA`, `NEA`, `STATUS`, `EAF`, `BETA`, `SE`, `P` and `N`. `EA` and `NEA` are the effect and non-effect alleles. `STATUS` is GWASLab's [seven-digit status code](https://cloufield.github.io/gwaslab/StatusCode/): the first two digits record genome build, followed by one digit each for identifier checking, coordinate checking, allele standardisation, reference alignment, and palindromic-variant/indel handling. A `9` means the corresponding check was not performed. The validator also rejects duplicate headers, empty tables and rows with inconsistent field counts.

The sidecar records `summary_statistics_id`, trait type and declared prevalence metadata, build, ancestry, origin, source format/method/release/basename and SHA-256, internal producer identity when applicable, canonical filename/SHA-256/columns/variant count, transformation and serialization, harmonisation metadata, and the optional non-secret access constraint. An already-gzipped canonical candidate is copied byte-for-byte; an uncompressed candidate is gzip-compressed deterministically.

All build-specific GWASLab reference parameters default to unset because no compact bundled reference is scientifically adequate. With no references, raw output is still standardised for names, columns and allele roles. Supplying `--gwaslab_reference_fasta_grch37` or `--gwaslab_reference_fasta_grch38` enables reference-allele checks and flips; the corresponding `--gwaslab_rsid_vcf_*` enables rsID assignment, and `--gwaslab_strand_vcf_*` enables palindromic-strand inference. The declared genome build chooses the resource set per summary.

GWASLab drops variants that fail its sanity checks and duplicated variants. Its log is not published because timestamps and container-local paths make it non-reproducible provenance noise.

## Heritability

Individual-level heritability output is method-first and then analysis under `individual/`. These native estimator tables are not harmonised.

### GCTA GREML and GREML-LDMS

<details markdown="1">
<summary>Output files</summary>

[GCTA](https://yanglab.westlake.edu.cn/software/gcta/) `.hsq` output reports genetic and residual variance components, `V(G)/Vp` on the observed scale, its standard error, log likelihood, likelihood-ratio test, p-value and `n`, the sample count after genotype, phenotype and covariate intersection. For a binary row with `population_prevalence`, GCTA also reports the native liability-scale `V(G)/Vp_L` line.

- `heritability/individual/gcta_greml/<analysis_id>/`
  - `<analysis_id>.gcta_greml.hsq`: Native GCTA GREML variance-component estimates.
- `heritability/individual/gcta_greml_ldms/<analysis_id>/`
  - `<analysis_id>.gcta_greml_ldms.hsq`: Native GCTA GREML-LDMS variance-component estimates.

</details>

GREML uses one dense matrix. GREML-LDMS partitions variants by LD score and MAF. Its method-option defaults are a `200` kb LD-score region, `4` LD strata and MAF boundaries `[0,0.01,0.05,0.2,0.5]`; changing any of them produces a distinct reusable matrix family.

### LDAK estimators

<details markdown="1">
<summary>Output files</summary>

[LDAK](https://dougspeed.com/ldak/) REML reports fitted kinship components, likelihood statistics and sample counts. Haseman-Elston regression provides the selected trait's regression estimate. PCGC is a binary-trait route and requires `population_prevalence` for ascertainment-aware estimation.

- `heritability/individual/ldak_reml/<analysis_id>/`
  - `<analysis_id>.ldak_reml.reml`: Native LDAK REML estimates.
  - `<analysis_id>.ldak_reml.reml.liab`: Optional liability-scale REML estimates for a binary analysis with `population_prevalence`.
- `heritability/individual/ldak_he/<analysis_id>/`
  - `<analysis_id>.ldak_he.he`: Native LDAK Haseman-Elston regression estimates.
- `heritability/individual/ldak_pcgc/<analysis_id>/`
  - `<analysis_id>.ldak_pcgc.pcgc`: Native LDAK PCGC regression estimates.
  - `<analysis_id>.ldak_pcgc.pcgc.marginal`: Optional marginal PCGC estimates emitted by LDAK.

</details>

The LDAK kinship model defaults to `human_default` with `power: -0.25`. Set `model: custom` before supplying another power. `weights_policy: equal` is the default and explicitly ignores weights; `default` retains LDAK's native policy; `provided` requires the staged `weights` resource. `relatedness_filter` defaults to `false`. When HE or PCGC has covariates, the pipeline first adjusts the kinship matrix on the same analysis subset and covariates, then passes those covariates to the estimator so phenotype residualisation and matrix projection remain aligned.

### Summary-level LDSC H2 and RG

<details markdown="1">
<summary>Output files</summary>

[LDSC](https://github.com/CBIIT/ldsc) consumes each canonical summary through one content-addressed HapMap3 munging step. Unary H2 and ordered pairwise RG requests reuse that munged result when the summary identity, adapter contract and HapMap3 bytes are identical. Munged summaries are workflow intermediates and are not published.

- `requests/ldsc_h2/<request_id>/`
  - `native.observed.log`: Complete native observed-scale LDSC H2 log.
  - `native.liability.log`: Optional native liability-scale log for a binary summary declaring both sample and population prevalence.
  - `diagnostics.tsv`: Completion classification, native scales, regression-SNP count, LDSC diagnostics and unmodified warnings.
  - `provenance.json`: Summary identity, request and reference-bundle metadata, accepted native arguments, munging keys/diagnostics, native estimates, LDSC source revision/container, warnings, classification and artifact inventory.
- `heritability/ldsc_h2/<request_id>/heritability.tsv`: One normalized H2 row per native scale.
- `requests/ldsc_rg/<request_id>/`: The same observed/optional-liability native logs, diagnostics and provenance for the ordered summary pair.
- `heritability/ldsc_rg/<request_id>/heritability.tsv`: Ordered left/right marginal H2 rows with endpoint-specific scale attribution.
- `genetic_covariance/ldsc_rg/<request_id>/genetic_covariance.tsv`: Genetic covariance and standard error with ordered scale attribution.
- `genetic_correlation/ldsc_rg/<request_id>/genetic_correlation.tsv`: Observed-invocation genetic correlation, standard error, z score and p value.

</details>

Normalized results preserve the native values and classify successful completion as `estimable`, `estimable_with_warning` or `completed_nonestimable`; native warnings and boundary violations never cause clipping or method selection. A malformed or incomplete mandatory native log fails the request. The pipeline presents every requested method and does not rank or combine them. Observed-scale LDSC is always retained. Liability-scale output is added only when all binary endpoints in the request declare both prevalence values. In a mixed RG invocation, a quantitative endpoint remains `observed`, the binary endpoint is `liability`, and the ordered covariance scale is written as `observed_x_liability` or `liability_x_observed`. The provenance retains LDSC's native labels separately from these scientifically attributed normalized scales.

## Pairwise GCTA bivariate REML

<details markdown="1">
<summary>Output files</summary>

The dense route uses one explicit all-variant GCTA matrix. The relationship's deterministic LDMS request owns one LD-by-MAF-stratified MGRM family. Neither route inherits matrix settings from an endpoint's unary analysis, and scientifically identical unary and pair LDMS settings reuse one matrix family. Native results remain request-addressed, while three lightweight TSV views expose each native estimand family without selecting, aggregating or ranking a result.

- `requests/gcta_bivariate_reml/<request_id>/`
  - `native.hsq`: Complete native GCTA bivariate REML variance-component result.
  - `native.log`: Native command, version, convergence and sample-overlap log.
  - `diagnostics.tsv`: Full-union endpoint counts, native common/non-missing counts, convergence, residual covariance when estimated, whether that component was retained or dropped by an explicit native option or GCTA's native overlap rule, warnings and the Q43 completion classification.
  - `provenance.json`: Ordered endpoint identities, cohort, matrix kind/key/native basename/settings, effective prevalence, accepted native arguments, all parsed native component values, tool version, warnings, classification and artifact inventory.
- `heritability/gcta_bivariate_reml/<request_id>/heritability.tsv`: Left and right `V(G)/Vp` estimates on every scale the native result emitted.
- `genetic_covariance/gcta_bivariate_reml/<request_id>/genetic_covariance.tsv`: Native observed-scale `C(G)_tr12` estimate and standard error.
- `genetic_correlation/gcta_bivariate_reml/<request_id>/genetic_correlation.tsv`: Native ordered `rG` estimate and standard error.
- `requests/gcta_bivariate_reml_ldms/<request_id>/`: The same native, diagnostics and provenance artifact set for REML-LDMS.
- `heritability/gcta_bivariate_reml_ldms/<request_id>/heritability.tsv`: Ordered trait-specific `V(Gk)/Vp` rows for every native LDMS component `Gk` and scale emitted by GCTA.
- `genetic_covariance/gcta_bivariate_reml_ldms/<request_id>/genetic_covariance.tsv`: Native observed-scale `C(Gk)_tr12` estimate and standard error for every LDMS component.
- `genetic_correlation/gcta_bivariate_reml_ldms/<request_id>/genetic_correlation.tsv`: Native ordered `rGk` estimate and standard error for every LDMS component.

</details>

Successful native completion is classified as `estimable`, `estimable_with_warning` or `completed_nonestimable`. An explicit native nonconvergence, corrupt or incomplete mandatory output, or execution error is the fourth state, `failed`; it fails the request and the run rather than publishing a misleading normalized result. A warning or out-of-range native estimate is retained rather than clipped or discarded. The pipeline does not compare methods or choose a best result. For binary endpoints, a declared `population_prevalence` is passed only through GCTA's endpoint-aware `--reml-bivar-prevalence` interface; ordinary unary `--prevalence` is never used on this route. Liability-scale heritability rows appear only when GCTA itself emits the corresponding native `_L` component.

## LDAK summary-statistics heritability and correlation

<details markdown="1">
<summary>Output files</summary>

SumHer and SumCors are request-addressed thin wrappers over native LDAK 6.3. Every native result is retained; the normalized TSVs are additional views and never rank, aggregate or select a preferred method.

- `requests/ldak_sumher/<request_id>/`
  - `native.hers`, `native.cats`, `native.share`, `native.enrich`, `native.extra`, `native.cross`, `native.taus`: Native SumHer estimates and category results.
  - `native.labels`, `native.progress`, `native.overlap`, `native.log`: Native labels, progress, overlap diagnostics and captured execution log.
  - `native.hers.liab`, `native.cats.liab`, `native.factor`: Optional native liability-scale artifacts, present only when LDAK receives a complete binary-trait prevalence/ascertainment pair.
  - `diagnostics.tsv`, `provenance.json`: Completion classification, warnings, exact request/reference identities, effective native arguments, adapter evidence, native log likelihoods, tool/runtime identity and native-artifact inventory. The wrapper does not invent a universal parameter count or derive AIC when the native output does not report the model-specific parameter count.
- `heritability/ldak_sumher/<request_id>/heritability.tsv`: Observed-scale and, when emitted by LDAK, liability-scale SumHer estimates.
- `requests/ldak_sumcors/<request_id>/`
  - `native.cors`, `native.cors.full`, `native.labels`, `native.progress`, `native.overlap`, `native.log`: Complete native SumCors result family.
  - `native.cors.liab`: Optional native liability-scale pair result when both ordered binary endpoints have complete prevalence/ascertainment declarations.
  - `diagnostics.tsv`, `provenance.json`: Ordered endpoint identities, reference ownership, effective native arguments, adapter evidence, classification, warnings and native-artifact inventory.
- `heritability/ldak_sumcors/<request_id>/heritability.tsv`: Ordered endpoint heritability estimates emitted by SumCors.
- `genetic_covariance/ldak_sumcors/<request_id>/genetic_covariance.tsv`: Native ordered coheritability estimate and uncertainty.
- `genetic_correlation/ldak_sumcors/<request_id>/genetic_correlation.tsv`: Native ordered genetic-correlation estimate and uncertainty.

</details>

Successful native completion uses the same `estimable`, `estimable_with_warning` and `completed_nonestimable` vocabulary as the GCTA bivariate views. Non-estimable compact or real datasets remain visible with native evidence and `NA` normalized estimates. Malformed or incomplete mandatory native results fail the request rather than being reinterpreted as a scientific result.

## Quality control and optional prepared data

### Relatedness matrices

<details markdown="1">
<summary>Output files</summary>

The pipeline builds each distinct [GCTA](https://yanglab.westlake.edu.cn/software/gcta/) or [LDAK](https://dougspeed.com/ldak/) relatedness matrix once and reuses it across compatible analyses. Matrices are unpublished by default because they are large intermediates. `<key>` is a content-derived digest over cohort identity and every construction setting that changes a matrix. Key-addressing is necessary because one cohort may need several matrices, while one matrix may belong to several analyses. Per-partition construction files and logs are never published.

- `quality_control/relatedness_matrices/<key>/` (with `--save_relatedness_matrices`)
  - `*.grm.bin`, `*.grm.N.bin`, `*.grm.id`: Dense GCTA matrix bundle.
  - `*.grm.sp`, `*.grm.id`: Sparse GCTA fastGWA matrix bundle.
  - `*.mgrm`, `*.grm.bin`, `*.grm.N.bin`, `*.grm.id`: GCTA GREML-LDMS manifest and stratified matrix bundles.
  - `*.grm.bin`, `*.grm.id`, `*.grm.details`, `*.grm.adjust`: LDAK matrix bundle.

</details>

Relatedness matrices are the only current `quality_control/` publication family; validation failures are reported before execution and do not create a published validation report.

### Prepared genotypes

<details markdown="1">
<summary>Output files</summary>

[PLINK 2](https://www.cog-genomics.org/plink/2.0/) converts PLINK 1 or VCF inputs once per cohort into the pipeline's canonical bundle. Only bundles the pipeline actually built are published. A cohort supplied as PLINK 2 is passed through without a process invocation, so it does not appear under `genotypes/`: the researcher's original `pgen`/`psam`/`pvar` already is the canonical bundle. Consequently, a four-cohort run with one PLINK 2 input can legitimately publish three prepared bundles; this does not mean a cohort was dropped.

- `genotypes/<cohort_id>/` (with `--save_prepared_genotypes`)
  - `<cohort_id>.pgen`: PLINK 2 genotype data converted by the pipeline.
  - `<cohort_id>.psam`: PLINK 2 sample information converted by the pipeline.
  - `<cohort_id>.pvar`: PLINK 2 variant information converted by the pipeline.

</details>

### Normalised phenotypes and covariates

<details markdown="1">
<summary>Output files</summary>

- `phenotypes/<analysis_id>/` (with `--save_normalised_phenotypes`)
  - `<analysis_id>.pheno`: Headered normalised phenotype file.
  - `<analysis_id>.qcovar`: Headered quantitative covariates, when supplied.
  - `<analysis_id>.catcovar`: Headered categorical covariates, when supplied.
  - `<analysis_id>.covar`: Headered merged covariates, when either covariate input was supplied.

</details>

These files show the exact recoding consumed by downstream tools and are useful for auditing case/control normalisation. They are unpublished by default because they are derived intermediates. Headerless tool-specific serialisations, the LDAK matrix-adjustment serialisation and the normalisation log are never published.

## MultiQC

<details markdown="1">
<summary>Output files</summary>

[MultiQC](https://multiqc.info/) combines the validated Analysis plan, workflow parameters, route-aware Methods Description and collected software versions into one report. The Analysis plan has one row per `analysis_id` and records the joined cohort, trait, trait type, genome build, ancestry provenance and requested association and heritability methods; the Methods Description additionally cites selected pairwise routes. It describes requested routes, not their completion or scientific results. Custom Manhattan and QQ plots and estimator-result panels are outside the current reporting scope.

- `multiqc/`
  - `multiqc_report.html`: Standalone HTML run report.
  - `multiqc_data/`: Parsed report data and the inputs used to build the report.
  - `multiqc_plots/`: Optional exported static plots.

</details>

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

[Nextflow](https://www.nextflow.io/docs/latest/reports.html) produces the execution reports, and nf-core records the launch parameters and versions used by the run. The version YAML contains the tools that actually executed; a tool absent because its route was not selected has no entry. MultiQC's own version is deliberately not fed into this upstream versions file because MultiQC consumes that file, which would create a circular channel dependency.

- `pipeline_info/`
  - `execution_report_<timestamp>.html`: Nextflow task and resource report.
  - `execution_timeline_<timestamp>.html`: Nextflow execution timeline.
  - `execution_trace_<timestamp>.txt`: Nextflow task trace.
  - `pipeline_dag_<timestamp>.html`: Nextflow workflow graph.
  - `pipeline_report.html`, `pipeline_report.txt`: Optional completion reports written when email reporting is requested.
  - `params_<timestamp>.json`: Complete run parameters.
  - `nf_core_gwas_software_mqc_versions.yml`: Software versions collected from executed processes.

</details>
