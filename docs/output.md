# nf-core/gwas: Output

## Introduction

This document describes the files that nf-core/gwas publishes beneath `--outdir`. This release prepares each declared cohort and each declared cohort-trait analysis and records how it did so; it runs no association or heritability method yet. Run-level provenance is collected in MultiQC and `pipeline_info/`.

Intermediates are unpublished by default. The optional directories described below appear only when their corresponding save control is enabled.

### Naming and attribution

Every published directory is attributed to one manifest identity:

- `genotypes/<cohort_id>/` is attributed to a `cohort_id` from the cohort manifest and carries that cohort's view record.
- `phenotypes/<analysis_id>/` is attributed to an `analysis_id` from the analysis manifest and carries the prepared inputs of that one cohort-trait analysis.

Use the following provenance chain for any published file:

1. Read `cohort_id` or `analysis_id` from its parent directory. Find that row in the retained manifest; an analysis row names its cohort through `cohort_id`.
2. Read the cohort's `genotypes/<cohort_id>/<cohort_id>.genotype_view.json` to learn which bytes the cohort's identity was taken from and which representation was on disk.
3. Read tool versions from `pipeline_info/nf_core_gwas_software_mqc_versions.yml`. The pipeline version and complete run parameters are recorded by the `pipeline_info/` reports and `params_<timestamp>.json`.

Together, the retained cohort and analysis manifests, the optional method-options document, the view records and the `pipeline_info/` artifacts identify the cohort, trait, genome build, pipeline version and producing tool version of every published file. Preserve them with an archived result.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and publishes:

- [Prepared genotypes](#prepared-genotypes)
- [Prepared phenotypes and covariates](#prepared-phenotypes-and-covariates)
- [MultiQC](#multiqc)
- [Pipeline information](#pipeline-information)

## Prepared genotypes

<details markdown="1">
<summary>Output files</summary>

Preparation preserves the representation a cohort was supplied in. A PLINK 1 or PLINK 2 cohort is used exactly as given, so nothing is built for it and nothing appears under `genotypes/` for it but its view record. [PLINK 2](https://www.cog-genomics.org/plink/2.0/) derives one shared PLINK 1 BED/BIM/FAM view of a PGEN cohort when — and only when — a selected method reads PLINK 1; no method is registered in this release, so no cohort is projected. A projection is a work-directory intermediate and is never published: what is published instead is the record of how it was made.

- `genotypes/<cohort_id>/`
  - `<cohort_id>.genotype_view.json`: The cohort's genotype view record, always written. `source` names the manifest's own format, the view identity (a `declared` token or a `sha256` digest of the supplied bytes), and the per-member basenames, sizes and digests where the pipeline computed them. `native_view` names the representation actually on disk and its reuse key. `plink1` names the PLINK 1 view: for a PLINK 1 cohort it is the supplied bundle itself and discards nothing, for a projected one it names its parent view, the full projection policy and the `information_loss` that policy accepts — `dosage` and `phase`. It is `null` for a cohort with no PLINK 1 consumer, which in this release is every cohort.

</details>

A cohort that declares `genotype_view_id` in its manifest row shows `source.identity` as `{"mode": "declared", "value": "<genotype_view_id>"}` and no member digests.

The view record is what will explain every other published key of that cohort once analysis routes are wired: an artifact built from a cohort's genotypes is shared by every cohort whose genotypes are byte-identical in the same representation, and its published directory is named by the artifact key rather than by a cohort — so when two `cohort_id`s name the same files, one directory serves both, and the two view records show why.

## Prepared phenotypes and covariates

<details markdown="1">
<summary>Output files</summary>

- `phenotypes/<analysis_id>/` (with `--save_normalised_phenotypes`)
  - `<analysis_id>.pheno`: Headered prepared phenotype file.
  - `<analysis_id>.qcovar`: Headered quantitative covariates, when supplied.
  - `<analysis_id>.catcovar`: Headered categorical covariates, when supplied.
  - `<analysis_id>.covar`: Headered merged covariates, when either covariate input was supplied.

</details>

These files show the exact representation that downstream analysis tools will consume: two identifier columns followed by the selected trait, with binary source values recoded to `0` and `1` and missing or unmatched values written as `NA`. They are unpublished by default because they are derived intermediates. Headerless tool-specific serialisations and the dummy-coded adjustment serialisation are never published. The public save parameter remains `--save_normalised_phenotypes`.

## MultiQC

<details markdown="1">
<summary>Output files</summary>

[MultiQC](https://multiqc.info/) combines the validated Analysis plan, workflow parameters, Methods Description and collected software versions into one report. The Analysis plan has one row per `analysis_id` and records the joined cohort, trait, trait type, genome build, ancestry provenance and requested association and heritability methods, which are `-` in this release. It describes the run plan, not the completion or scientific results of any step.

- `multiqc/`
  - `multiqc_report.html`: Standalone HTML run report.
  - `multiqc_data/`: Parsed report data and the inputs used to build the report.
  - `multiqc_plots/`: Optional exported static plots.

</details>

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

[Nextflow](https://www.nextflow.io/docs/latest/reports.html) produces the execution reports, and nf-core records the launch parameters and versions used by the run. The version YAML contains the tools that actually executed; a tool absent because no task of it ran has no entry. MultiQC's own version is deliberately not fed into this upstream versions file because MultiQC consumes that file, which would create a circular channel dependency.

- `pipeline_info/`
  - `execution_report_<timestamp>.html`: Nextflow task and resource report.
  - `execution_timeline_<timestamp>.html`: Nextflow execution timeline.
  - `execution_trace_<timestamp>.txt`: Nextflow task trace.
  - `pipeline_dag_<timestamp>.html`: Nextflow workflow graph.
  - `pipeline_report.html`, `pipeline_report.txt`: Optional completion reports written when email reporting is requested.
  - `params_<timestamp>.json`: Complete run parameters.
  - `nf_core_gwas_software_mqc_versions.yml`: Software versions collected from executed processes.

</details>
