# nf-core/gwas: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/gwas, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [#95](https://github.com/nf-core/gwas/pull/95) - Add the `regenie/runl1` module, completing the REGENIE Step 1 module set (`runl0`, `splitl0`, `runl1`, `step1`, `step2`).
- [#93](https://github.com/nf-core/gwas/pull/93) - Added LDAK local modules (addgrms, adjustgrm, calcgenotypeerrort2, calcinflation, calckins, createthinweights, filterrelatedness, he, kvikstep1, kvikstep2, pcgc, reml, thinpredictors) under `modules/local/ldak/`.

### `Fixed`

- [#96](https://github.com/nf-core/gwas/pull/96) - Bump the `nft-utils` nf-test plugin to `0.0.9` for nf-test 0.9.4 compatibility (fixes `sanitizeOutput()` `MissingMethodException` in module tests after the 4.0.2 template merge).

### `Dependencies`

- Updated the `plink/gwas` (`2d5c9c0`) and `plink/vcf` (`6d46786`) modules to the latest nf-core/modules version, preserving the pipeline's local patches (the `plink/vcf` `--pheno`/`--make-bed` phenotype input and the `plink/gwas` split `assoc`/`qassoc` outputs plus resource-label overrides).

### `Deprecated`
