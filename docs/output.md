# nf-core/gwas: Output

## Introduction

This document describes the files that nf-core/gwas publishes beneath `--outdir`. Native association, heritability and declared pairwise results are retained. Every internal association result and every external summary source passes through GWASLab, whose table and log are published directly. Run-level provenance is collected in MultiQC and `pipeline_info/`. Native results are retained unchanged.

Intermediates are unpublished by default. The optional directories described below appear only when their corresponding save control is enabled.

### Naming and attribution

The shared result prefix grammar is:

```text
<analysis_id>.<method>[.<shard>]
```

`<analysis_id>` is copied from the analysis manifest and identifies one cohort-trait analysis unit. `<method>` is one of the method-selector tokens documented in [Usage](usage.md#relational-manifest-input). A producing tool can add a native result suffix after that prefix.

Summary results use a separate first-class identity. Pipeline-generated association summaries use `<analysis_id>--<association_method>`; external summaries use the declared `summary_statistics_id`. Both publish as `<summary_statistics_id>.gwaslab.tsv.gz` with `<summary_statistics_id>.gwaslab.log` in the same identity-addressed directory.

Use the following provenance chain for any result:

1. For an analysis result, read `<analysis_id>` and `<method>` from its parent directories and filename. Find that analysis row, follow its `cohort_id`, and inspect its method options.
2. For a summary result, read `summary_statistics_id` from its directory and filename. Resolve its declared internal producer or external source from the retained summary-statistics manifest.
3. Map the method to its producing tool using the table below.
4. Read tool versions from `pipeline_info/nf_core_gwas_software_mqc_versions.yml`. The pipeline version and complete run parameters are recorded by the `pipeline_info/` reports and `params_<timestamp>.json`.

Pairwise outputs instead use the deterministic request ID `<method>--<relationship_id>`. Find `relationship_id` in `--relationship_manifest`, follow its ordered left and right endpoint IDs, and inspect the native result and log under `requests/<method>/<request_id>/`. The pipeline does not add a normalised estimand table or a diagnostics table, and never rewrites a native file. It publishes each native result unchanged.

| Method token                                                                                                                                      | Producing tool |
| ------------------------------------------------------------------------------------------------------------------------------------------------- | -------------- |
| `regenie`                                                                                                                                         | REGENIE        |
| `gcta_fastgwa`, `gcta_greml`, `gcta_greml_ldms`, `gcta_bivariate_reml`, `gcta_bivariate_reml_ldms`, `gcta_bivariate_he`, `gcta_bivariate_he_ldms` | GCTA           |
| `ldak_kvik`, `ldak_reml`, `ldak_he`, `ldak_pcgc`, `ldak_fast_he`, `ldak_fast_pcgc`                                                                | LDAK 6         |
| `mph_reml`, `mph_reml_ldms`, `mph_bivariate_reml`, `mph_bivariate_reml_ldms`                                                                      | MPH 0.55.1     |
| `ldak_sumher`, `ldak_sumcors`                                                                                                                     | LDAK 6.3       |
| `ldsc_h2`, `ldsc_rg`                                                                                                                              | LDSC           |

Together, the result prefix, retained cohort and analysis manifests, optional method-options document, and `pipeline_info/` artifacts identify the analysis, cohort, trait, genome build, method, scientific settings, pipeline version and producing tool version. Preserve them with an archived result.

Analysis attribution applies under `association/` and `heritability/individual/`; summary attribution applies under `summary_statistics/` and summary-level request outputs. Genotype views and relatedness matrices are deliberately shared artifacts rather than trait-method results: `genotypes/` is attributed to `cohort_id` and carries each cohort's view record, while `quality_control/relatedness_matrices/` and `quality_control/ldms_component_plans/` are attributed to their own reuse keys and may each serve several analysis rows — including rows on different `cohort_id`s whose genotypes are byte-identical.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and publishes:

- [Association](#association)
  - [REGENIE](#regenie)
  - [GCTA fastGWA-MLM](#gcta-fastgwa-mlm)
  - [LDAK-KVIK](#ldak-kvik)
- [Summary statistics](#summary-statistics)
- [Heritability](#heritability)
  - [GCTA GREML and GREML-LDMS](#gcta-greml-and-greml-ldms)
  - [LDAK estimators](#ldak-estimators)
  - [LDAK direct-genotype estimators](#ldak-direct-genotype-estimators)
  - [MPH REML and REML-LDMS](#mph-reml-and-reml-ldms)
  - [Summary-level LDSC H2 and RG](#summary-level-ldsc-h2-and-rg)
- [Pairwise GCTA bivariate REML and HEreg](#pairwise-gcta-bivariate-reml-and-hereg)
- [Pairwise MPH bivariate REML](#pairwise-mph-bivariate-reml)
- [LDAK summary-statistics heritability and correlation](#ldak-summary-statistics-heritability-and-correlation)
- [Quality control and optional prepared data](#quality-control-and-optional-prepared-data)
- [MultiQC](#multiqc)
- [Pipeline information](#pipeline-information)

## Association

Native association output is published by method and then analysis. The native tables are not rewritten; deterministic filename tokens added by tools are removed where necessary so every published name follows the shared prefix grammar.

### REGENIE

<details markdown="1">
<summary>Output files</summary>

[REGENIE](https://rgcgithub.github.io/regenie/) fits one whole-genome prediction model per analysis in Step 1 and tests variants in Step 2. Analyses retain separate fits even when their inputs are identical. Standard Step 1 is the default because it is the simplest execution path; `--regenie_step1_mode chunked` and `--regenie_step1_jobs` use split-L0/run-L0/run-L1 when a large cohort needs work divided into smaller jobs. `--regenie_lowmem` defaults to `true` so temporary prediction blocks stay in the task work directory rather than memory.

- `association/regenie/<analysis_id>/`
  - `<analysis_id>.regenie.gz`: Native, space-delimited REGENIE Step 2 association result.

</details>

For binary traits, the per-analysis `regenie.firth`, `regenie.firth_approx` and `regenie.firth_p_threshold` method options default to approximate Firth correction below `0.01`. `regenie.min_mac` defaults to `null`, retaining REGENIE's own versioned minimum-MAC policy unless an analysis explicitly overrides it.

The published Step 2 file is native output with the fixed `_PHENO` token removed from its name. Step 1 predictions, chunk-planning files, temporary low-memory predictions and logs remain in the work directory. Resumed Step 1 reuse therefore requires preserving the Nextflow cache and work outputs.

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

[LDAK-KVIK](https://dougspeed.com/ldak-kvik/) fits a Step 1 prediction model from the PLINK 1 compatibility bundle prepared once per cohort and tests the full bundle in Step 2. `ldak.kvik_step1_subset` defaults to `all`; `thin_common` requests deterministic thinning and `provided` requires the stageable `ldak.predictor_extract` resource. One `thin_common` predictor artifact is built per compatible prepared genotype view and effective thinning contract, then shared by every Step 1 fit on that view, so phenotype and covariate differences do not fragment thinning reuse. Each analysis declares its own Step 1 fit using the resolved predictor artifact: two analyses that resolve to the same predictor artifact still fit one model each. The native `.assoc` file is published unchanged. The internal three-file Step 1 bundle is consumed directly by Step 2; predictor lists, effects, thinning progress and logs remain in the work directory. Resumed Step 1 reuse therefore requires preserving the Nextflow cache and work outputs.

- `association/ldak_kvik/<analysis_id>/`
  - `<analysis_id>.ldak_kvik.step2.assoc`: Native LDAK-KVIK Step 2 association table.

</details>

## Summary statistics

<details markdown="1">
<summary>Output files</summary>

[GWASLab](https://cloufield.github.io/gwaslab/) standardises every pipeline-generated association result and every external source; native association results remain available under `association/`. External rows declare one explicit GWASLab input format, including `gwaslab` for a pre-harmonised GWASLab table. No source bypasses the process, and the external source itself is not republished.

- `summary_statistics/<summary_statistics_id>/`
  - `<summary_statistics_id>.gwaslab.tsv.gz`: Gzip-compressed, tab-delimited GWASLab-standard summary-statistics table.
  - `<summary_statistics_id>.gwaslab.log`: GWASLab harmonisation log for the same source.

</details>

The table uses GWASLab's standard field names, including `SNPID`, `CHR`, `POS`, `EA`, `NEA`, `STATUS`, `EAF`, `BETA`, `SE`, `P` and `N` where the declared source format supplies or derives them. `EA` and `NEA` are the effect and non-effect alleles. `STATUS` is GWASLab's [seven-digit status code](https://cloufield.github.io/gwaslab/StatusCode/). Downstream program adapters consume this GWASLab artifact directly; there is no second canonical serializer, post-GWASLab schema validator or custom checksum/provenance artifact.

All build-specific GWASLab reference parameters default to unset because no compact bundled reference is scientifically adequate. With no references, raw output is still standardised for names, columns and allele roles. Supplying `--gwaslab_reference_fasta_grch37` or `--gwaslab_reference_fasta_grch38` enables reference-allele checks and flips; the corresponding `--gwaslab_rsid_vcf_*` enables rsID assignment, and `--gwaslab_strand_vcf_*` enables palindromic-strand inference. The declared genome build chooses the resource set per summary.

GWASLab drops variants that fail its sanity checks and duplicated variants. The emitted log is published beside the table as the native record of harmonisation.

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
  - `<analysis_id>.ldak_he.he.liab`: Optional liability-scale HE estimates for a binary analysis with `population_prevalence`.
  - `<analysis_id>.ldak_he.factor`: Native observed-to-liability conversion factors when prevalence is supplied.
- `heritability/individual/ldak_pcgc/<analysis_id>/`
  - `<analysis_id>.ldak_pcgc.pcgc`: Native LDAK PCGC regression estimates.
  - `<analysis_id>.ldak_pcgc.pcgc.marginal`: Optional marginal PCGC estimates emitted by LDAK.

</details>

The LDAK kinship model defaults to `human_default` with `power: -0.25`. Set `model: custom` before supplying another power. `weights_policy: equal` is the default and applies no predictor weights; `provided` requires the staged `weights` resource passed to native `--weights`. `relatedness_filter` defaults to `false`. When HE or PCGC has covariates, the pipeline first adjusts the kinship matrix on the same analysis subset and covariates, then passes those covariates to the estimator so phenotype residualisation and matrix projection remain aligned.

### LDAK direct-genotype estimators

<details markdown="1">
<summary>Output files</summary>

[LDAK](https://dougspeed.com/ldak/) fast Haseman-Elston and fast PCGC estimate the same quantities as `ldak_he` and `ldak_pcgc` but read the cohort's PLINK 1 genotypes directly, approximating the kinship trace terms with random vectors (the randomised approach introduced by [RHE-mc](https://doi.org/10.1038/s41467-020-17576-9)). No relatedness matrix is built, requested or published for these routes, so they add nothing under `quality_control/`. Covariates are projected out natively from the complete covariate files preparation writes.

- `heritability/individual/ldak_fast_he/<analysis_id>/`
  - `<analysis_id>.ldak_fast_he.fasthe`: Native randomised Haseman-Elston estimates: a key-value header (`Num_Kinships`, `Num_Top_Predictors`, `Num_Covariates`, `Coeffsfile`, `Covar_Heritability`, `Total_Samples`, `With_Phenotypes`, `Null_Likelihood`, `Alt_Likelihood`, `LRT_Stat`, `LRT_P`) followed by a `Component Heritability SE Size Mega_Intensity SE` table with `Her_K1`, `Her_Top` and `Her_All` rows.
  - `<analysis_id>.ldak_fast_he.fasthe.liab`: Optional liability-scale FASTHE estimates for a binary analysis with `population_prevalence`.
  - `<analysis_id>.ldak_fast_he.factor`: Native observed-to-liability conversion factors when prevalence is supplied.
  - `<analysis_id>.ldak_fast_he.cats`, `.share`, `.enrich`, `.cross`: Native per-category heritability, share, enrichment and cross-product results.
  - `<analysis_id>.ldak_fast_he.coeff`: Native covariate effects (`Component Effect SE P`), one `Covariate_n` row per fitted column plus the intercept.
  - `<analysis_id>.ldak_fast_he.log`: Native log recording the complete argument list, the effective random-vector count, the jackknife block count, the seed when one was supplied, the sample and predictor counts, and any weights-coverage warning.
- `heritability/individual/ldak_fast_pcgc/<analysis_id>/`
  - `<analysis_id>.ldak_fast_pcgc.fastpcgc`: Native randomised PCGC estimates in the same layout, on the liability scale.
  - `<analysis_id>.ldak_fast_pcgc.fastpcgc.marginal`: Native marginal PCGC estimate, unconditional on the covariates and therefore different from the primary result whenever covariates were fitted.
  - `<analysis_id>.ldak_fast_pcgc.cats`, `.share`, `.enrich`, `.cross`, `.coeff`, `.log`: As above, with `.coeff` reporting `Component Log_Odds SE P`.

</details>

The per-predictor (`.ind.hers`), per-block (`.jackests`), per-random-vector (`.repetitions`), label, progress and combined-covariate files stay in the work directory: they are large or purely diagnostic. Three things about these results are worth stating plainly. `ldak.fast_num_blocks` sets the number of **predictor** jackknife blocks — LDAK partitions predictors, not samples — so it moves the standard error and not the point estimate, and a block-jackknife standard error is not comparable across block settings nor to the exact `ldak_he`/`ldak_pcgc` standard error. An unseeded run cannot be reproduced: LDAK records no seed in its log unless one was supplied, so set `ldak.fast_seed` for anything you intend to publish. And with `weights_policy: provided`, inspect the `.log` for `contains weights for only`, which means LDAK gave weight zero to every predictor the weights file omitted.

### MPH REML and REML-LDMS

<details markdown="1">
<summary>Output files</summary>

[MPH](https://jiang18.github.io/mph/) fits a REML variance-component model by MINQUE, estimating the trace terms with random vectors, over relationship matrices it builds itself. `mph_reml` fits one component over the autosomal variant universe; `mph_reml_ldms` fits one component per LD-by-MAF stratum of the shared component plan. Both add the residual component `err`, both are quantitative-only, and both publish MPH's own files unchanged.

- `heritability/individual/mph_reml/<analysis_id>/` and `heritability/individual/mph_reml_ldms/<analysis_id>/`
  - `<analysis_id>.<method>.mq.vc.csv`: Primary result. `trait_x`, `trait_y`, `vc_name`, `m`, `var`, `seV`, `pve`, `seP`, `enrichment`, `seE`, one row per component in matrix-list order and `err` last, followed by the enrichment and variance-component sampling-covariance blocks. `vc_name` is the staged matrix prefix verbatim and carries the stratum key for a stratified fit; `m` is that component's post-quality-control variant count and is `NA` on the `err` row. Variances are unconstrained, so a negative one is native. The appended covariance columns repeat their names, so a header-keyed reader will silently collapse them.
  - `<analysis_id>.<method>.mq.blue.csv`: Best linear unbiased estimates of the fitted covariates: `trait`, `covar`, `blue`, `se`, `pval`, then their sampling covariance in columns named `<trait>.<covar>`. `covar` names the columns MPH actually fitted.
  - `<analysis_id>.<method>.mq.iter.csv`: Solver trace, one row per completed iteration: `iter`, `num_traits`, `sample_size`, `num_GRMs`, `logLL`, `dLLpred`, `dogleg_Newton`. `sample_size` is the analysis set the fit used.
  - `<analysis_id>.<method>.log`: Complete MPH standard output: the echoed option block including the effective thread count, the `Non-missing analysis set contains N individuals` line, the random-vector line and the per-iteration trust-region trace. MPH writes no summary of its estimates anywhere and prints its two quality failures — non-convergence and a rank-deficient covariate design — only here, as `Warning:` lines at exit 0.

</details>

`<analysis_id>.<method>.mq.py.csv`, MPH's fixed-effect-adjusted phenotype, is deliberately not published: it is one row per individual and grows with the cohort, so it stays in the task work directory for the same reason LDAK's per-individual heritabilities do. The matrix list the fit is addressed by stays there too.

Two things about these results are worth stating plainly. A seed makes an MPH fit repeatable, not stable: measured on a 200-sample cohort over seeds 1 to 8, the proportion of variance explained spanned `0.052` at MPH's default of 100 random vectors and `0.114` at 50, while the reported standard error moved only between `0.124` and `0.129` and therefore contains none of that Monte-Carlo component — so set `mph.random_vectors` explicitly for anything you intend to publish. And an estimate is reproducible only against the same seed, random-vector count, thread count and memory mode: the last two move it in the sixth to seventh significant digit.

### Summary-level LDSC H2 and RG

<details markdown="1">
<summary>Output files</summary>

[LDSC](https://github.com/CBIIT/ldsc) consumes each GWASLab summary through one content-addressed HapMap3 munging step. Unary H2 and ordered pairwise RG requests reuse that munged result when the summary identity, adapter contract and HapMap3 bytes are identical. Munged summaries are workflow intermediates and are not published.

- `requests/ldsc_h2/<request_id>/`
  - `native.log`: Complete native LDSC H2 log, on the liability scale for a binary summary declaring both sample and population prevalence, otherwise on the observed scale.
- `requests/ldsc_rg/<request_id>/`
  - `native.log`: Complete native LDSC RG log, with liability conversion when at least one endpoint is binary and every binary endpoint declares both prevalence values.

</details>

The pipeline presents every requested LDSC invocation without converting its log into a common heritability, covariance or correlation family. Each request runs once. Prevalence flags select native liability conversion only when at least one endpoint is binary and all binary endpoints in the request declare both prevalence values; otherwise the run uses the observed scale. Quantitative endpoints use LDSC's native `nan` placeholder in a mixed RG invocation.

## Pairwise GCTA bivariate REML and HEreg

<details markdown="1">
<summary>Output files</summary>

The dense route uses one explicit all-variant GCTA matrix. The relationship's deterministic LDMS request owns one LD-by-MAF-stratified MGRM family. Neither route inherits matrix settings from an endpoint's unary analysis, and scientifically identical unary and pair LDMS settings reuse one matrix family regardless of which method token requested it. `gcta_bivariate_he` and `gcta_bivariate_he_ldms` run GCTA's `--HEreg-bivar` Haseman-Elston cross-product (HE-CP) estimator on that same dense matrix or MGRM family. HE-CP only makes the fitting stage cheaper than REML; it still requires the same full dense (or LDMS-stratified) matrix construction, storage and I/O, so it is published as a deterministic moment reference and sensitivity analysis, not a matrix-free or more scalable route. Each method retains its separate native result contract.

- `requests/gcta_bivariate_reml/<request_id>/`
  - `native.hsq`: Complete native GCTA bivariate REML variance-component result.
  - `native.log`: Native command, version, convergence and sample-overlap log.
- `requests/gcta_bivariate_reml_ldms/<request_id>/`
  - `native.hsq`: Complete native GCTA bivariate REML-LDMS result: one `V(Gk)_tr1`, `V(Gk)_tr2` and `C(Gk)_tr12` row per stratum, the residual rows, the phenotypic variances, and one genetic correlation `rGk` per stratum. GCTA does not report a genome-wide total for a multi-component bivariate REML fit.
  - `native.log`: Native command, version, convergence and sample-overlap log. It also carries the complete sampling variance/covariance matrix of the variance-component estimates, which GCTA writes only here, and — when the request declared a `population_prevalence` — the liability-scale `V(Gk)/Vp_tr<k>_L` rows, which the bivariate `.hsq` does not carry although the unary `.hsq` does.
- `requests/gcta_bivariate_he/<request_id>/`
  - `native.HEreg`: Complete native GCTA `--HEreg-bivar` dense result: `Intercept_tr1`, `Intercept_tr2`, `Intercept_tr12`, `V(G)/Vp_tr1`, `V(G)/Vp_tr2`, `C(G)/Vp_tr12`, `rG`, `N_tr1` and `N_tr2` rows, each with `Estimate`, `SE_OLS`, `SE_Jackknife`, `P_OLS` and `P_Jackknife` columns.
  - `native.log`: Native command and version, plus the complete jackknife sampling variance/covariance matrix of the estimates; GCTA writes that matrix only to the log, never to `.HEreg`.
- `requests/gcta_bivariate_he_ldms/<request_id>/`
  - `native.HEreg`: Native GCTA HEreg-LDMS table, including its per-component and native total rows.
  - `native.log`: Native command, version and jackknife sampling variance/covariance matrix.

</details>

The pipeline preserves GCTA's native values and does not compare methods or choose a preferred result. For binary endpoints, a declared `population_prevalence` is passed only through GCTA's endpoint-aware `--reml-bivar-prevalence` interface; ordinary unary `--prevalence` is never used on this route. `gcta_bivariate_he` and `gcta_bivariate_he_ldms` accept quantitative pairs only and do not accept pair covariates. HEreg fits the declared left-trait-by-right-trait orientation on the lower triangle of the GRM, so reversed duplicate relationships remain invalid.

`gcta_bivariate_reml_ldms` publishes only the per-stratum correlations in GCTA's native result. GCTA does not report a genome-wide genetic correlation for its multi-component bivariate REML model, and the pipeline does not derive one.

## Pairwise MPH bivariate REML

<details markdown="1">
<summary>Output files</summary>

[MPH](https://jiang18.github.io/mph/) fits the two declared endpoints jointly by REML, estimating the trace terms with random vectors, over relationship matrices it builds itself. `mph_bivariate_reml` fits one component over the autosomal variant universe; `mph_bivariate_reml_ldms` fits one component per LD-by-MAF stratum of the shared component plan, whose settings the relationship's deterministic LDMS request owns. Both add the residual component `err`, both accept quantitative pairs only, and both publish MPH's own files unchanged.

- `requests/mph_bivariate_reml/<request_id>/` and `requests/mph_bivariate_reml_ldms/<request_id>/`
  - `native.mq.vc.csv`: Primary variance-component result. `trait_x`, `trait_y`, `vc_name`, `m`, `var`, `seV`, `pve`, `seP`, `enrichment`, `seE`, one row per component in matrix-list order and `err` last, for each trait and for the cross pair, followed by the enrichment and variance-component sampling-covariance blocks. `vc_name` is the staged matrix prefix verbatim and carries the stratum key for a stratified fit; `m` is that component's post-quality-control variant count and is `NA` on the `err` row. Variances and covariances are unconstrained, so a negative one is native. The appended covariance columns repeat their names, so a header-keyed reader will silently collapse them.
  - `native.mq.cor.csv`: Native correlations: `vc_name`, `trait_x`, `trait_y`, `cor` and `se`, with one row per fitted component, a residual row `err` and the genome-wide row `G`.
  - `native.mq.blue.csv`: Best linear unbiased estimates of the fitted covariates, one block per trait: `trait`, `covar`, `blue`, `se`, `pval`, then their sampling covariance in columns named `<trait>.<covar>`. `covar` names the columns MPH actually fitted.
  - `native.mq.iter.csv`: Solver trace, one row per completed iteration: `iter`, `num_traits`, `sample_size`, `num_GRMs`, `logLL`, `dLLpred`, `dogleg_Newton`. `sample_size` is the analysis set the fit used.
  - `native.log`: Complete MPH standard output: the echoed option block including the effective thread count, the `Non-missing analysis set contains N individuals` line, the random-vector line and the per-iteration trust-region trace. MPH writes no summary of its estimates anywhere and prints its quality failures only here, as `Warning:` lines at exit 0.

</details>

The genome-wide total is MPH's own number. MPH prints the genome-wide genetic correlation and its standard error itself, as the synthetic `G` row of `native.mq.cor.csv`, for the one-component and the multi-component fit alike; for `mph_bivariate_reml` that row is the single component's row. The pipeline publishes what MPH wrote and derives nothing, so these two routes have no derived estimand table at all. GCTA's `gcta_bivariate_reml_ldms` route publishes only native per-stratum correlations, while MPH fits the complete-case intersection of the pair and GCTA's bivariate REML fits an unbalanced design over the same two endpoints.

MPH does not always form a correlation for every row. Its variance components are unconstrained, so a component can carry a negative fitted genetic variance for one of the two traits, and that component's correlation then comes out `-nan` with a `nan` standard error, at exit 0, with every other row of the same file complete. This happens on the shipped compact fixture. The native file retains MPH's missing-value representation unchanged.

Read the native pair by trait name rather than by column position. The declared left/right order is written into MPH's `--trait_names` in that order, and MPH then labels its output pair in reverse — `--trait_names A,B` produces rows labelled `trait_x = B`, `trait_y = A` — consistently, in `native.mq.cor.csv` and in the cross block of `native.mq.vc.csv`.

MPH fits the complete-case intersection: both traits and every named covariate must be observed for an individual to enter the fit. GCTA's bivariate REML keeps every individual over the same partially overlapping pair and fits an unbalanced design instead, so the two published results are not matched models and their difference is not attributable to the estimator alone. The native log and solver trace report the analysis-set size.

The sample is one of two things that separate the two pair families. The other is the constraint: MPH's REML is unconstrained and has no option to constrain, while the pipeline's GCTA bivariate routes run GCTA's constrained default, since their rendered arguments are the endpoint-aware prevalence and the request's own `native_args` and nothing more. On a pair where that constraint binds — a low-signal pair will do it — a published MPH result and a published GCTA result are two different estimators rather than two implementations of one, whatever their samples. `--reml-no-constrain` on the GCTA pair request is the setting that matches them, and it is the researcher's to add.

`native.mq.py.csv`, MPH's fixed-effect-adjusted phenotype, is deliberately not published: it is one projected phenotype per individual and grows with the cohort, so it stays in the task work directory for the same reason the unary route's does. The matrix list the fit is addressed by stays there too.

An MPH pair estimate is reproducible against the same seed, random-vector count, thread count and memory mode, and not otherwise. MPH's unset seed is the fixed integer `0` rather than entropy, so an unseeded run repeats; but the seed value moves the point estimate and not merely its standard error, and the thread count the executor chooses and `save_memory` move it in the sixth to seventh significant digit. The native log records the effective settings. Raise `random_vectors` on the pair request for anything you intend to publish, and treat both routes as a recommended REML candidate rather than an established default until parity with the GCTA pair routes has been demonstrated.

## LDAK summary-statistics heritability and correlation

<details markdown="1">
<summary>Output files</summary>

SumHer and SumCors are request-addressed wrappers over native LDAK 6.3. Every native result is retained directly; the pipeline does not create normalized estimand, diagnostics or provenance views.

- `requests/ldak_sumher/<request_id>/`
  - `native.hers`, `native.cats`, `native.share`, `native.enrich`, `native.extra`, `native.cross`, `native.taus`: Native SumHer estimates and category results.
  - `native.labels`, `native.progress`, `native.overlap`, `native.log`: Native labels, progress, overlap diagnostics and captured execution log.
  - `native.hers.liab`, `native.cats.liab`, `native.factor`: Optional native liability-scale artifacts, present only when LDAK receives a complete binary-trait prevalence/ascertainment pair.
- `requests/ldak_sumcors/<request_id>/`
  - `native.cors`, `native.cors.full`, `native.labels`, `native.progress`, `native.overlap`, `native.log`: Complete native SumCors result family.
  - `native.cors.liab`: Optional native liability-scale pair result when both ordered binary endpoints have complete prevalence/ascertainment declarations.

</details>

The native files retain LDAK's own result structure, warnings and missing-value representation. They are not parsed into a pipeline-wide estimand vocabulary.

## Quality control and optional prepared data

### Relatedness matrices

<details markdown="1">
<summary>Output files</summary>

The pipeline builds each reusable [GCTA](https://yanglab.westlake.edu.cn/software/gcta/), [LDAK](https://dougspeed.com/ldak/) or [MPH](https://jiang18.github.io/mph/) base artifact once and derives each requested child artifact once. Compatible GREML and fastGWA routes share a dense GCTA base, while the fastGWA cutoff identifies only its sparse child. Filtered and unrestricted LDAK routes share a kinship base, while unrelated-sample subsetting identifies only its child. LDAK HE and PCGC likewise share a covariate-adjusted child when their selected parent, effective sample subset, covariate content and native options match.

Matrices are unpublished by default because they are large intermediates. `<key>` is a content-derived artifact identity over immutable input or parent identity and effective scientific settings. It never contains focal analysis, estimator or publication state. Key-addressing is necessary because one cohort may need several artifacts, while one artifact may serve unary, pairwise and association consumers. Each base or child is published exactly once by its own key; per-partition construction files and logs are never published.

- `quality_control/relatedness_matrices/<key>/` (with `--save_relatedness_matrices`)
  - `*.grm.bin`, `*.grm.N.bin`, `*.grm.id`: Dense GCTA matrix bundle.
  - `*.grm.sp`, `*.grm.id`: Sparse GCTA fastGWA child bundle.
  - `*.grm.bin`, `*.grm.N.bin`, `*.grm.id`: GCTA GREML-LDMS stratified matrix bundles in their declared non-empty LD-by-MAF order. Each native consumer writes its small MGRM control list inside its own task; the control list is not an independently published artifact.
  - `*.grm.bin`, `*.grm.id`, `*.grm.details`, `*.grm.adjust`: Base or unrelated-sample child LDAK kinship bundle.
  - `*.grm.bin`, `*.grm.id`, `*.grm.details`, `*.grm.adjust`, `*.grm.root`: Covariate-adjusted LDAK child bundle.
  - `*.grm.bin`, `*.grm.iid`: One-component or LD-by-MAF-stratified MPH matrix bundle, the stratified family in its declared non-empty order. There is no `*.grm.N.bin`.

</details>

An MPH bundle and a GCTA bundle are never interchangeable, and renaming one into the other is not a conversion. MPH writes an 8-byte header holding the sample count and the sum of the SNP weights, an unnormalised row-major upper triangle, and an IID-only companion; GCTA writes a headerless normalised row-major lower triangle with `.grm.N.bin` beside it and a FID/IID companion. Neither tool refuses the other's layout usefully: a GCTA bundle relabelled for MPH runs to exit 0 and produces nothing, and an MPH bundle relabelled for GCTA is read without complaint — `gcta --pca` on one returned eigenvalues of 11262.9, 7958.5 and 7701.4 against the true 2.78, 2.67 and 2.59, at exit 0. Build the matrix each route asks for instead.

Relatedness matrices and LD-by-MAF component plans are the two current `quality_control/` publication families; validation failures are reported before execution and do not create a published validation report. LDAK fast Haseman-Elston and fast PCGC build no matrix and therefore publish nothing here.

### LDMS component plans

<details markdown="1">
<summary>Output files</summary>

An LD-by-MAF component plan is the LD scores of one cohort genotype view and the ordered, disjoint SNP groups derived from them. It is built once by GCTA's LD-score pass and the pipeline's stratifier, and every stratified matrix family that declares the same settings reads it, so `gcta_greml_ldms` and `mph_reml_ldms` on one row partition one identical variant set. `<plan_key>` digests that genotype view and the three plan settings only, never the tool or the method token that asked for it, which is why the plan is published as its own family rather than inside either matrix family.

- `quality_control/ldms_component_plans/<plan_key>/` (with `--save_relatedness_matrices`)
  - `*_gcta_ld.score.ld`: Native GCTA per-variant LD scores for the plan's declared region size.
  - `*.strata.tsv`: Ordered stratum manifest, LD-major and MAF-minor with empty strata omitted: `model_key`, `stratum_key`, `ld_lower`, `ld_upper`, `maf_lower`, `maf_upper`, `predictor_count`, `group_filename`.
  - `*_snp_group_<stratum_key>.txt`: The variant IDs of one non-empty stratum, one per line.

</details>

The manifest's row order is the component order of every matrix family and every result built from the plan, and `stratum_key` appears verbatim in the `vc_name` of an MPH stratified result, so a component can be traced from the published estimate back to its variant list by name rather than by row position. `predictor_count` is what the plan declares; the count a tool retained after its own quality control is reported by that tool's result.

### Prepared genotypes

<details markdown="1">
<summary>Output files</summary>

Preparation preserves the representation a cohort was supplied in. A PLINK 1 or PLINK 2 cohort is used exactly as given, so nothing is built for it and nothing appears under `genotypes/` for it but its view record. [PLINK 2](https://www.cog-genomics.org/plink/2.0/) derives one shared PLINK 1 BED/BIM/FAM view of a PGEN cohort when — and only when — a selected method reads PLINK 1. That projection is a work-directory intermediate and is never published: what is published instead is the record of how it was made.

- `genotypes/<cohort_id>/`
  - `<cohort_id>.genotype_view.json`: The cohort's genotype view record, always written. `source` names the manifest's own format, the view identity (a `declared` token or a `sha256` digest of the supplied bytes), and the per-member basenames, sizes and digests where the pipeline computed them. `native_view` names the representation actually on disk and its reuse key. `plink1` names the PLINK 1 view: for a PLINK 1 cohort it is the supplied bundle itself and discards nothing, for a projected one it names its parent view, the full projection policy and the `information_loss` that policy accepts — `dosage` and `phase`. It is `null` for a cohort with no PLINK 1 consumer.

</details>

A cohort that declares `genotype_view_id` in its manifest row, such as a [mounted VCZ cohort](usage.md#mounted-vcz-cohorts), shows `source.identity` as `{"mode": "declared", "value": "<genotype_view_id>"}` and no member digests.

The view record is what explains every other published key of that cohort. A relatedness matrix, a REGENIE Step 1 fit and an LDAK thin-common artifact are shared by every cohort whose genotypes are byte-identical in the same representation, and their published directories are named by the artifact key rather than by a cohort — so when two `cohort_id`s name the same files, one directory serves both, and the two view records show why.

### Prepared phenotypes and covariates

<details markdown="1">
<summary>Output files</summary>

- `phenotypes/<analysis_id>/` (with `--save_normalised_phenotypes`)
  - `<analysis_id>.pheno`: Headered prepared phenotype file.
  - `<analysis_id>.qcovar`: Headered quantitative covariates, when supplied.
  - `<analysis_id>.catcovar`: Headered categorical covariates, when supplied.
  - `<analysis_id>.covar`: Headered merged covariates, when either covariate input was supplied.

</details>

These files show the exact representation consumed by downstream tools. They are unpublished by default because they are derived intermediates. Headerless tool-specific serialisations and the LDAK matrix-adjustment serialisation are never published. The public save parameter remains `--save_normalised_phenotypes`.

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

### Meta-analysis results

A meta-analysis row publishes `<summary_statistics_id>.gwaslab.tsv.gz` below its ordinary
`summary_statistics/<summary_statistics_id>/` directory. This derived table retains fixed-effect variants
meeting `min_studies` and appends native fields from each explicitly selected model. `BETA`, `SE` and `P`
always describe the GWASLab inverse-variance fixed-effect result. Its native `EAF`, `N`, `Z`, `Q`, `DOF`,
`P_HET`, `I2` and `DIRECTION` fields are retained; `N_STUDIES` is `DOF + 1`.

| Selected model | Additional fields in the derived table                                                                                                                  |
| -------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `random`       | `BETA_RANDOM`, `SE_RANDOM`, `Z_RANDOM`, `P_RANDOM`                                                                                                      |
| `re2`          | `P_RE2`, `RE2_MEAN_COMPONENT`, `RE2_HET_COMPONENT`                                                                                                      |
| `mrmega`       | `MRMEGA_CHISQ_ASSOC`, `MRMEGA_DF_ASSOC`, `MRMEGA_P_ASSOC`; corresponding `CHISQ`, `DF`, `P` fields for `ANCESTRY_HET` and `RESIDUAL_HET`; `MRMEGA_LNBF` |

Optional model fields are joined to the retained fixed-effect rows by allele-aware variant identity. An absent
native model row leaves those optional fields missing. MR-MEGA's multi-degree-of-freedom association p-value
stays in `MRMEGA_P_ASSOC`, separate from the pooled fixed-effect estimate. Native numerical spelling, zero
p-values and `NA` statistics are preserved; no tail probabilities are recomputed.

Native products are also retained unchanged:

- `requests/fixed/<request_id>/`: `<request_id>.fixed.tsv.gz` and `<request_id>.gwaslab.log`.
- `requests/random/<request_id>/`: `<request_id>.random.tsv.gz` when selected.
- `requests/re2/<request_id>/`: `<request_id>.metasoft.txt` and `<request_id>.metasoft.log` when selected.
- `requests/mrmega/<request_id>/`: `<request_id>.result` and `<request_id>.log` when selected, including native
  coefficients and per-marker diagnostics.

The native files contain each program's full output before the derived-summary contribution threshold.
METASOFT's ordinary fixed and random fields remain inside its native file and are not separate published
model families. Aligned study matrices and MR-MEGA input lists remain in the Nextflow work directory.
