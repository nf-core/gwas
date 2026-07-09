# nf-core/gwas: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/gwas, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

- Updated the `plink/gwas` (`2d5c9c0`) and `plink/vcf` (`6d46786`) modules to the latest nf-core/modules version, preserving the pipeline's local patches (the `plink/vcf` `--pheno`/`--make-bed` phenotype input and the `plink/gwas` split `assoc`/`qassoc` outputs plus resource-label overrides).

### `Deprecated`
