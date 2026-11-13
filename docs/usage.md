# nf-core/gwas: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/gwas/usage](https://nf-co.re/gwas/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/gwas takes genotyped cohorts and the traits to analyse on them as two linked manifests. A cohort manifest owns genotype facts; an analysis manifest links each individual-level trait analysis to one cohort. This release validates the two manifests, prepares each cohort's genotypes and each analysis's phenotype and covariates once, and reports the run plan; the association and heritability routes that consume those prepared inputs are added in later releases, each registering its own method selector token.

> [!IMPORTANT]
> Genotypes must be prepared before you run the pipeline. The pipeline converts accepted genotype
> encodings into the formats required by its methods, but it does not perform genotype quality control.

## Relational manifest input

Every run supplies the linked `--cohort_manifest` and `--analysis_manifest` together.

The JSON schemas own column names, required fields, types, enumerations and file existence. Schema-optional columns may be omitted entirely or supplied with blank cells. Unexpected and repeated column names are rejected; quoted CSV headers and values are supported, and column order is not significant.

```bash
nextflow run nf-core/gwas \
    -r <VERSION> \
    -profile docker \
    --cohort_manifest ./cohorts.csv \
    --analysis_manifest ./analyses.csv \
    --outdir ./results
```

The cohort manifest owns reusable genotype facts. The analysis manifest owns one cohort–trait analysis and joins to exactly one cohort through `cohort_id`. Each `cohort_id` has one non-conflicting definition; each `analysis_id` is unique; and every analysis references a declared cohort. Cohort preparation is shared by every analysis on the cohort, while outputs retain `analysis_id`.

### Cohort manifest fields

| Column             | Required | Description                                                                                                                                                                                                                              |
| ------------------ | -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `cohort_id`        | Yes      | Unique, whitespace-free cohort identifier and analysis-manifest foreign key.                                                                                                                                                             |
| `genome_build`     | Yes      | `GRCh37` or `GRCh38`; recorded in provenance and selects build-specific resources once a route consumes them.                                                                                                                            |
| `ancestry`         | Yes      | Case-sensitive provenance label beginning with a letter or digit and containing only letters, digits, `_`, `.`, or `-`.                                                                                                                  |
| `pgen`             | By group | PLINK 2 genotype file; supply the complete `pgen`/`psam`/`pvar` group.                                                                                                                                                                   |
| `psam`             | By group | PLINK 2 sample file.                                                                                                                                                                                                                     |
| `pvar`             | By group | PLINK 2 variant file ending in `.pvar`. A Zstandard-compressed `.pvar.zst` is not accepted: genotype consumers that address a PLINK 2 fileset by a single prefix read only an uncompressed `<prefix>.pvar`.                              |
| `bed`              | By group | PLINK 1 genotype file; supply the complete `bed`/`bim`/`fam` group.                                                                                                                                                                      |
| `bim`              | By group | PLINK 1 variant file.                                                                                                                                                                                                                    |
| `fam`              | By group | PLINK 1 sample file.                                                                                                                                                                                                                     |
| `genotype_view_id` | No       | Immutable identity of the supplied genotype view: 16–64 lowercase hexadecimal characters, optionally prefixed by a backend namespace such as `biofuse:`. When absent, the pipeline digests the supplied bytes once per cohort in a task. |

Populate exactly one complete genotype representation on each row, and give the members of a multi-file group one shared basename stem — genotype consumers that address a fileset by a single prefix require it. VCF input is not accepted; convert a VCF to PGEN or BED with PLINK 2 before declaring the cohort.

The supplied representation is preserved. A PLINK 1 or PLINK 2 cohort is used exactly as given and nothing is converted for it. A PLINK 1 BED/BIM/FAM view is derived from a PGEN cohort once, and only when a selected method reads PLINK 1; no method is registered in this release, so no cohort is projected. That projection takes hard calls at an explicit `--hard-call-threshold 0.1` and drops dosage and phase; every cohort's `genotypes/<cohort_id>/<cohort_id>.genotype_view.json` records which of these applied to it.

### Analysis manifest fields

| Column                  | Required | Description                                                                                                                                                         |
| ----------------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `analysis_id`           | Yes      | Unique, whitespace-free result identifier.                                                                                                                          |
| `cohort_id`             | Yes      | A `cohort_id` declared in the cohort manifest.                                                                                                                      |
| `trait_id`              | Yes      | Whitespace-free trait identifier retained in provenance.                                                                                                            |
| `trait_type`            | Yes      | `quantitative` or `binary`; never inferred from values.                                                                                                             |
| `phenotype`             | Yes      | Existing headered phenotype file containing the selected trait.                                                                                                     |
| `phenotype_column`      | Yes      | Header name of the selected trait column.                                                                                                                           |
| `control_value`         | Binary   | Source value recoded to `0`; required for binary traits, forbidden for quantitative traits, and distinct from `case_value`.                                         |
| `case_value`            | Binary   | Source value recoded to `1`; required for binary traits and forbidden for quantitative traits.                                                                      |
| `quant_covariates`      | No       | Existing headered quantitative-covariate file.                                                                                                                      |
| `cat_covariates`        | No       | Existing headered categorical-covariate file.                                                                                                                       |
| `association_methods`   | No       | Comma-delimited association selector. No association method is registered in this release, so the cell must be empty; each route adds its token when it is wired.   |
| `heritability_methods`  | No       | Comma-delimited heritability selector. No heritability method is registered in this release, so the cell must be empty; each route adds its token when it is wired. |
| `population_prevalence` | By route | Number strictly between `0` and `1`; valid only for a binary analysis whose selected estimator consumes it. With no method registered, it is rejected on every row. |
| `sample_prevalence`     | No       | Optional sample case fraction strictly between `0` and `1` for binary traits; forbidden for quantitative traits.                                                    |

Tokens in a selector cell are comma-delimited without spaces and may appear only once; a token that no registered method claims is rejected by name, and the diagnostic lists the accepted values. A row with both selector cells empty declares a cohort-trait pairing to prepare, which is every row in this release.

### Examples

[`assets/examples/relational/cohort_manifest.csv`](../assets/examples/relational/cohort_manifest.csv) declares the public two-autosome PGEN cohort from the nf-core/test-datasets `gwas` branch. The [preparation-only analysis manifest](../assets/examples/relational/analysis_manifest_preparation_only.csv) declares one quantitative trait on it with both selector cells empty, and is what `-profile test` runs:

```bash
nextflow run nf-core/gwas \
    -r <VERSION> \
    -profile docker \
    --cohort_manifest assets/examples/relational/cohort_manifest.csv \
    --analysis_manifest assets/examples/relational/analysis_manifest_preparation_only.csv \
    --save_normalised_phenotypes \
    --outdir results
```

The [quantitative example](../assets/examples/relational/analysis_manifest_quantitative.csv) is the same row with an association selector populated. It shows the shape of a method-selecting row and is rejected by this release, which registers no method.

## Advanced method options

`--method_options` is optional. The document keeps its JSON root keyed by `analysis_id`, or places that same object under an `analyses` namespace; each analysis entry names a method family and the options it sets for it. No method family is registered in this release, so every family is rejected by name and an empty document is the only accepted one; each route adds its family, option names and defaults when it is wired. Unlisted analysis settings receive their defaults.

### Validation diagnostics

Structural failures name the manifest and invalid column. Cross-row preflight failures additionally name the CSV row and `cohort_id` or `analysis_id`: examples include an incomplete genotype group, a second genotype group, members that do not share one basename stem, conflicting duplicates, an undefined cohort, unknown or repeated method tokens, invalid binary coding and route-inapplicable prevalence. Method-options failures name the document, analysis and qualified option before task submission.

### Phenotype input preparation

`PREPARE_PHENOTYPE_INPUTS` prepares each selected phenotype in a shared two-identifier-plus-trait layout. Binary source values matching `control_value` and `case_value` become `0` and `1`; missing values and unmatched binary values become `NA`. Quantitative values must be numeric. Quantitative and categorical covariates remain distinct for tools with separate native interfaces, and a merged covariate file is written when either was supplied.

### Run-level defaults

Parameters divide by what they control. Scientific specification — what is estimated and under which model — is declared per analysis in the method-options document, not on the command line. CPU, memory and time requests are not parameters at all: change them with a custom process configuration passed via `-c`.

| Parameter or behaviour | Default and rationale                                                                                                           |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| Save controls          | Off. Intermediates stay out of the results directory unless explicitly requested, avoiding unexpectedly large published output. |

The one opt-in save control is:

- `--save_normalised_phenotypes`: publish headered prepared phenotype and covariate files under `phenotypes/<analysis_id>/`. The parameter name is retained for compatibility.

Prepared genotype derivatives have no save control: a cohort's supplied files are never copied, and a derived PLINK 1 projection is a work-directory intermediate described by the always-published view record.

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

Pin the pipeline revision with `-r`, retain the exact supplied manifests, method-options document and parameter file, and archive `pipeline_info/` with your results. Reusing the same revision, inputs and parameters also lets `-resume` recover cached tasks.

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

The pipeline assigns default CPU, memory and time requests by process label and retries selected resource-related failures once, with doubled requests. Cap every request at what your machine or queue offers with `process.resourceLimits` in a custom configuration file, and [customise individual processes](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) with `withName` selectors.

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

| What you see                                             | Meaning                                                                                      | First evidence to inspect                                                                                                |
| -------------------------------------------------------- | -------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| An error naming a manifest, CSV row, identity and field  | Preflight rejected input before task submission                                              | Read the complete grouped error report, fix every listed row and relaunch.                                               |
| `Method-options document ... analysis_id ... option ...` | The method-options document names an analysis, family or option the pipeline does not accept | Read the analysis, option and reason; then correct the document.                                                         |
| `Process ... terminated with an error`                   | A task was submitted and failed                                                              | Inspect `.nextflow.log` and the task's `.command.err`, `.command.out` and `.command.log` in the reported work directory. |

### Manifest validation failed

Read the reported CSV path, row number, entity ID, field name and reason from left to right. The validator reports all linked-manifest errors it can find in one launch, so correct every bullet before rerunning.

| Diagnostic fragment                                                                                                       | Supported action                                                                                                                                                          |
| ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Validation of samplesheet failed` and a named field                                                                      | Correct the schema violation. Required columns must be present; schema-optional columns may be omitted; unexpected columns are rejected. Column order is not significant. |
| `header row 1 repeats column ...`                                                                                         | Keep each column name once. Quoted header names are accepted.                                                                                                             |
| `no genotype group is populated`, `a second genotype group is populated` or `genotype group ... is only partly populated` | Populate exactly one complete representation: `pgen`/`psam`/`pvar` or `bed`/`bim`/`fam`.                                                                                  |
| `do not share one basename stem`                                                                                          | Rename the members of the group so that they share one prefix.                                                                                                            |
| `duplicate cohort_id` or `duplicate analysis_id`                                                                          | Keep one row for each identity. Conflicting duplicates also report the fields that differ.                                                                                |
| `undefined cohort_id`                                                                                                     | Make the analysis row's `cohort_id` match one cohort-manifest identity exactly.                                                                                           |
| `unknown method` or `listed more than once`                                                                               | Use only the accepted values the diagnostic lists, each at most once. This release lists none, so leave both selector cells empty.                                        |
| A binary `case_value` or `control_value` diagnostic                                                                       | Supply both distinct source codes for a binary trait. Remove both from a quantitative row. The pipeline does not infer `trait_type` from phenotype values.                |
| A `population_prevalence` diagnostic                                                                                      | Use it only for a binary analysis whose selected estimator consumes it. See [Analysis manifest fields](#analysis-manifest-fields).                                        |

### Method-options validation failed

The diagnostic names the document, analysis ID, fully qualified option and reason. Fix malformed JSON, name only analyses declared in the analysis manifest, and use only the `analyses` namespace described under [Advanced method options](#advanced-method-options). Operational settings such as process resources and tool threads belong in run or profile configuration, not `--method_options`.

### A process failed

1. Copy the process name and work-directory hash from the terminal or `.nextflow.log`.
2. Inspect `.command.err`, `.command.out` and `.command.log` in that task directory; `.command.sh` records the executed command.
3. Correct the reported input, configuration, container, filesystem or resource problem.
4. Relaunch the same command with `-resume`; do not delete `work/` or `.nextflow/` first.

Selected resource failures are retried with larger requests, subject to `process.resourceLimits`. See [Resource requests](#resource-requests) and the general [nf-core troubleshooting guide](https://nf-co.re/docs/running/troubleshooting).

### Expected results are missing

- Prepared phenotypes and covariates are unpublished by default; enable `--save_normalised_phenotypes` before expecting `phenotypes/`.
- Every cohort publishes its view record under `genotypes/<cohort_id>/`; the genotype files themselves are never copied.
- Published intermediates are retention outputs, not importable cross-run caches. Supported reuse requires retained Nextflow work and `-resume`.

### `-resume` did not reuse work

Nextflow reuses a task only when its inputs, pipeline code and relevant configuration still match a retained cache entry. Keep the original `work/` directory and `.nextflow/` cache, rerun from the same launch context where practical, and use `nextflow log` to find a run name for `-resume <run-name>`. A changed input file, manifest, method option, pipeline revision, process configuration or missing work directory can require recomputation. See [`-resume`](#-resume).

### Get help

Before requesting help, retain the pipeline version, exact launch command with secrets removed, `.nextflow.log`, the failing process name and relevant `.command.*` files, the reported manifest row or a minimal redacted reproducer, executor/profile and container runtime, and whether the original work directory remains available for `-resume`.

Use the [nf-core troubleshooting guide](https://nf-co.re/docs/running/troubleshooting), open a [GitHub issue](https://github.com/nf-core/gwas/issues) with a reproducible report, or ask in the nf-core [`#gwas` Slack channel](https://nfcore.slack.com/channels/gwas).
