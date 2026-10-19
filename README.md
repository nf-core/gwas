<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-gwas_logo_dark.png">
    <img alt="nf-core/gwas" src="docs/images/nf-core-gwas_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/gwas)
[![GitHub Actions CI Status](https://github.com/nf-core/gwas/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/gwas/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/gwas/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/gwas/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/gwas/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.6-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.3-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.3)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/gwas)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23gwas-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/gwas)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/gwas** is a bioinformatics pipeline for association, individual-level and summary-level heritability, declared pairwise genetic-correlation analysis, and common-variant meta-analysis. A cohort manifest owns genotype facts, an analysis manifest links traits and individual-level methods to those cohorts, a summary-statistics manifest declares external, pipeline-generated or meta-analysis-derived summary results and their unary methods, and an optional relationship manifest binds explicit analysis or summary endpoints. The pipeline reuses compatible work, retains native analytical results, and passes every internal and external summary-statistics source through one GWASLab harmonisation boundary.

Genotype quality control is not performed by the pipeline. Input genotypes must already have suitable samples, variants, alleles, coordinates, genome build and analysis filters.

![Overview of the nf-core/gwas analysis routes](docs/images/nf-core-gwas_metro_map.svg)

## Pipeline summary

1. Validate linked cohort, analysis, summary-statistics, relationship, reference and request declarations.
2. Prepare each distinct cohort once from PLINK 2 or PLINK 1 input, preserving the supplied representation.
3. Prepare the selected phenotype and optional quantitative and categorical covariates for the native tools.
4. Run selected association routes:
   - REGENIE
   - GCTA fastGWA-MLM
   - LDAK-KVIK
5. Assign every association result a stable `<analysis_id>--<association_method>` summary-statistics identity and standardise it with GWASLab while preserving the native result.
6. Load every external summary-statistics table through GWASLab using an explicit named format, including pre-harmonised tables declared with the `gwaslab` format.
7. Build reusable base relatedness artifacts once and derive sparse, subset or adjusted children for selected association and heritability routes:
   - GCTA GREML
   - GCTA GREML-LDMS
   - LDAK REML
   - LDAK Haseman-Elston regression
   - LDAK PCGC
   - MPH REML
   - MPH REML-LDMS

   The two LD- and MAF-stratified routes share one component plan: the LD scores of a cohort's genotype view and the ordered, disjoint SNP groups derived from them are built once per view and settings, and both the GCTA and the MPH stratified matrix families are built from that one plan. MPH builds its matrices in its own layout, which is not interchangeable with GCTA's.

   Not every heritability route needs one: LDAK fast Haseman-Elston regression and LDAK fast PCGC estimate directly from the cohort's genotypes and build no relatedness matrix at all.

8. Pool ordered base summaries with fixed effects and optional conventional random effects, METASOFT RE2 and MR-MEGA meta-regression; retain native model results and make the derived fixed-effect summary available to downstream summary methods.
9. Resolve declared unary and pair summary-statistics requests against explicit LDAK or LDSC reference bundles, including LDAK SumHer heritability and SumCors genetic correlation.
10. Run declared same-cohort pairs with dense or LDMS GCTA bivariate REML or HEreg, or with dense or LDMS MPH bivariate REML, and publish each native result and log directly. The ordered two-trait table and the relationship-owned covariates are prepared once per relationship and shared by every selected pair estimator, so several estimators over one relationship still define one endpoint sample set.
11. Run declared LDSC H2 and ordered RG requests, reusing content-identical munging and publishing the observed- and available liability-scale native logs directly.
12. Collect run and software provenance with MultiQC and Nextflow reports.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

Prepare a cohort manifest and an analysis manifest linked by `cohort_id`. Advanced per-analysis
scientific settings and stageable resources may be supplied in an optional method-options JSON
document.

```csv title="cohorts.csv"
cohort_id,genome_build,ancestry,pgen,psam,pvar,bed,bim,fam
my_cohort,GRCh37,EUR,/data/my_cohort.pgen,/data/my_cohort.psam,/data/my_cohort.pvar,,,
```

```csv title="analyses.csv"
analysis_id,cohort_id,trait_id,trait_type,phenotype,phenotype_column,control_value,case_value,quant_covariates,cat_covariates,association_methods,heritability_methods,population_prevalence,sample_prevalence
height,my_cohort,height,quantitative,/data/phenotypes.tsv,height,,,,,regenie,,,
```

Runnable minimal and heterogeneous examples are available under [`assets/examples/relational/`](assets/examples/relational/).

To ingest external summaries or select unary summary methods for a pipeline-generated association result, add the summary-statistics manifest. To request pairwise analysis, add the optional relationship manifest; no pair is inferred. Individual-level pairs accept two distinct analysis IDs from the same cohort and select any of `gcta_bivariate_reml`, `gcta_bivariate_reml_ldms`, `gcta_bivariate_he`, `gcta_bivariate_he_ldms`, `mph_bivariate_reml` and `mph_bivariate_reml_ldms`. Summary requests select one named LDAK or LDSC bundle through `--method_options`, with the staged resource roles declared once in `--reference_catalog`.

Then run:

```bash
nextflow run nf-core/gwas \
    -r <VERSION> \
    -profile docker \
    --cohort_manifest cohorts.csv \
    --analysis_manifest analyses.csv \
    --relationship_manifest relationships.csv \
    --outdir results
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/gwas/usage) and the [parameter documentation](https://nf-co.re/gwas/parameters).

## Pipeline output

Native association results are published under `association/<method>/<analysis_id>/`. Every internal or external summary result is published once under `summary_statistics/<summary_statistics_id>/` as `<summary_statistics_id>.gwaslab.tsv.gz` with its GWASLab log. Individual-level heritability estimates remain under `heritability/individual/<method>/<analysis_id>/`. Declared GCTA, MPH, LDAK and LDSC requests publish their native result families under `requests/<method>/<request_id>/`, and the pipeline assembles no common estimand family. Every native result is retained unchanged. A per-result invocation record is added for the MPH routes, whose runtime does not record its own invocation. Munged LDSC summaries and other adapters are not republished; prepared genotypes, prepared phenotypes and relatedness matrices remain unpublished unless their save controls are enabled. REGENIE and LDAK-KVIK Step 1 predictions stay in Nextflow work and are consumed directly by Step 2.

For exact filenames, output layout and optional prepared data, see the [output documentation](https://nf-co.re/gwas/output).

## Credits

nf-core/gwas was originally written by Chris Wyatt, Fernando Duarte, Maxime Laurent.

We thank the following people for their extensive assistance in the development of this pipeline:

- [lyh970817](https://github.com/lyh970817)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#gwas` channel](https://nfcore.slack.com/channels/gwas) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add the pipeline citation when a Zenodo DOI is available, then update the badge above. -->
<!-- If you use nf-core/gwas for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
