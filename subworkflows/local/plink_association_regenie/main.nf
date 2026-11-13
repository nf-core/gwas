// Fit explicitly declared REGENIE producers and attach each complete prediction bundle to its test requests.
// Fit and test genotypes are independent; test results retain their analysis metadata and native shard files.

include { PLINK_FIT_REGENIE } from '../plink_fit_regenie/main'
include { REGENIE_STEP2     } from '../../../modules/local/regenie/step2/main'

workflow PLINK_ASSOCIATION_REGENIE {
    take:
    ch_fits // channel: [ val(fit_id), val(meta), path(genotype), path(variants), path(samples), path(pheno), path(covar), val(bsize), val(mode), val(jobs) ]; covar and standard-mode jobs may be []
    ch_predictions // channel: [ val(fit_id), val(meta), path(predictions), path(loco) ]; supplied producers, disjoint from ch_fits
    ch_tests // channel: [ val(fit_id), val(meta), path(genotype), path(variants), path(samples), path(pheno), path(covar), val(bsize) ]; one bundle per test event

    main:
    // The nested fitter joins its inputs on meta.id. Its private metadata uses the explicit producer key;
    // the side channel restores the caller's complete metadata on the public prediction artifact.
    ch_fit_metadata = ch_fits.map { fit_id, meta, _genotype, _variants, _samples, _pheno, _covar, _bsize, _mode, _jobs -> tuple(fit_id, meta) }
    ch_fit_inputs = ch_fits.multiMap { fit_id, meta, genotype, variants, samples, pheno, covar, bsize, mode, jobs ->
        def fit_meta = meta + [id: fit_id]
        genotypes: tuple(fit_meta, genotype, variants, samples)
        pheno: tuple(fit_meta, pheno)
        covar: tuple(fit_meta, covar)
        bsize: tuple(fit_meta, bsize)
        mode: tuple(fit_meta, mode)
        jobs: tuple(fit_meta, jobs)
    }

    PLINK_FIT_REGENIE(ch_fit_inputs.genotypes, ch_fit_inputs.pheno, ch_fit_inputs.covar, ch_fit_inputs.bsize, ch_fit_inputs.mode, ch_fit_inputs.jobs)

    ch_produced = PLINK_FIT_REGENIE.out.predictions
        .map { meta, predictions -> tuple(meta.id, predictions) }
        .join(PLINK_FIT_REGENIE.out.loco.map { meta, loco -> tuple(meta.id, loco) }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_fit_metadata, failOnDuplicate: true, failOnMismatch: true)
        .map { fit_id, predictions, loco, meta -> tuple(fit_id, meta, predictions, loco) }
    ch_bundles = ch_produced.mix(ch_predictions)

    // A producer may serve many analyses and chromosome bundles. Combine on F only, never on A or a
    // genotype basename; this stages every companion LOCO file for every native Step 2 invocation.
    ch_step2 = ch_tests
        .combine(ch_bundles, by: 0)
        .multiMap { _fit_id, meta, genotype, variants, samples, pheno, covar, bsize, fit_meta, predictions, loco ->
            genotypes: tuple(meta, genotype, variants, samples)
            predictions: tuple(fit_meta, predictions, loco)
            pheno: tuple(meta, pheno)
            covar: tuple(meta, covar)
            bsize: bsize
        }
    REGENIE_STEP2(ch_step2.genotypes, ch_step2.predictions, ch_step2.pheno, ch_step2.covar, ch_step2.bsize)

    // The nested fitter reports under the private producer key; restore the caller's metadata as on the
    // prediction artifact so a fit stage stays diagnosable from this boundary.
    ch_fit_logs = PLINK_FIT_REGENIE.out.logs
        .map { meta, log -> tuple(meta.id, log) }
        .combine(ch_fit_metadata, by: 0)
        .map { fit_id, log, meta -> tuple(fit_id, meta, log) }

    emit:
    predictions = ch_bundles // channel: [ val(fit_id), val(meta), path(predictions), path(loco) ]; once per explicit producer
    results     = REGENIE_STEP2.out.results // channel: [ val(meta), path(results) ]; native files per test event, no chromosome merge
    logs        = REGENIE_STEP2.out.log // channel: [ val(meta), path(log) ]; test-event logs
    fit_logs    = ch_fit_logs // channel: [ val(fit_id), val(meta), path(log) ]; native logs of every fit stage this call ran
}
