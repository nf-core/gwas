# LDAK SumCors test fixtures

These compact text fixtures are deterministically derived from the official LDAK [Test Datasets](https://dougspeed.com/test-datasets/) human PLINK dataset and its `quant.pheno` and `quant2.pheno` phenotypes.

The first 512 predictors whose allele pair is not A/T, T/A, C/G, or G/C were selected in `human.bim`. The checksum-pinned LDAK 6.3 compatibility image was then run with:

```console
ldak --linear trait1 --bfile human --pheno quant.pheno --extract predictors.txt --max-threads 2
ldak --linear trait2 --bfile human --pheno quant2.pheno --extract predictors.txt --max-threads 2
ldak --calc-tagging test --bfile human --extract predictors.txt --power -0.25 --max-threads 2
```

Only `trait1.summaries`, `trait2.summaries`, and `test.tagging`, the native mandatory inputs exercised by this module, are retained. They should move to `nf-core/test-datasets` before an upstream module submission.
