# ![nf-core/gwas](docs/images/nf-core-gwas_logo.png)

**This is a new skeleton pipeline**.

[![GitHub Actions CI Status](https://github.com/nf-core/gwas/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/gwas/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/gwas/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/gwas/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/gwas.svg)](https://hub.docker.com/r/nfcore/gwas)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23gwas-4A154B?logo=slack)](https://nfcore.slack.com/channels/gwas)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Podman`](https://podman.io/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/gwas -profile test,<docker/singularity/podman/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```bash
    nextflow run nf-core/gwas -profile <docker/singularity/podman/conda/institute> --input samplesheet.csv --genome GRCh37
    ```

See [usage docs](https://nf-co.re/gwas/usage) for all of the available options when running the pipeline.

## Documentation

The nf-core/gwas pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/gwas/usage) and [output](https://nf-co.re/gwas/output).

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

<!-- From Raquel Genomics imputaiton pipeline-->


A tool for imputation of genotype array datasets from dbGaP. The Genotype Imputation Pipeline consists of the following steps:

0. Identify input genome build version automatically
1. Lift the input to build GRCh37 (hg19)
2. Quality control 1: LD-based fix of strand flips, fix strand swaps, filter variants by missingness
3. Split samples by ancestry
4. Quality control 2: filter samples by missingness, filter variants by HWE
5. Phase
6. Impute


## Dependencies
The pipeline was tested in garibaldi using the following required software and packages:

- R v3.5.1
- vcftools v0.1.14
- PLINK v1.9
- PLINK v2.00a3LM 64-bit Intel
- samtools v1.9
- GenotypeHarmonizer v1.4.20
- ADMIXTURE
- Eagle v2.4
- Minimac4
- liftOver

## How to run


### Step 0: Check genome build and select chain file


Where:
- myinput is the full path to the input genotype array dataset in either vcf or vcf.gz format
- myoutput is the full path to save the output of this step
- gz (gz=yes or gz=no) is whether the input file is either vcf or vcf.gz format

The output file will have the sufix *.BuildChecked

### Step 1: Lifeover input genotype array to GRCh37 build

### Step 2: LD-based fix of strand flips, fix strand swaps and mismatching alleles, and initial quality control (90% missingnes per variant)

### Step 3: Estimate ancestry and split samples by ancestry

### Step 4: 2nd quality control

### Step 5: Phasing

### Step 6: Imputation and post-imputation quality control



## Credits

nf-core/gwas was originally written by Abhinav Sharma, Raquel Dias, Alain Coletta.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#gwas` channel](https://nfcore.slack.com/channels/gwas) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/gwas for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
