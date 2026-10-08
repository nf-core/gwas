// Fit declared KVIK producers and attach complete predictions to explicitly referenced testing requests.
include { LDAK_THINCOMMON } from '../../../modules/local/ldak/thincommon/main'
include { LDAK_KVIKSTEP1  } from '../../../modules/local/ldak/kvikstep1/main'
include { LDAK_KVIKSTEP2  } from '../../../modules/local/ldak/kvikstep2/main'

workflow PLINK_ASSOCIATION_LDAK_KVIK {
    take:
    ch_fits // channel: [ val(fit_id), val(meta), path(bed), path(bim), path(fam), path(phenotype), path(qcovar), path(covar), path(extract), val(thin_common) ]
    ch_predictions // channel: [ val(fit_id), val(meta), path(root), path(details), path(prs) ], supplied producers bypass fitting
    ch_tests // channel: [ val(fit_id), val(meta), path(bed), path(bim), path(fam), path(phenotype), path(qcovar), path(covar), path(keep) ]

    main:
    // The native producer uses F privately; the original caller metadata is restored on every fit output.
    def ch_fit_metadata = ch_fits.map { fit_id, meta, _bed, _bim, _fam, _phenotype, _qcovar, _covar, _extract, _thin_common -> [fit_id, meta] }
    def ch_fit_records = ch_fits.map { fit_id, meta, bed, bim, fam, phenotype, qcovar, covar, extract, thin_common ->
        [fit_id, meta + [id: fit_id], bed, bim, fam, phenotype, qcovar, covar, extract, thin_common]
    }
    def ch_fit_branches = ch_fit_records.branch { _fit_id, _meta, _bed, _bim, _fam, _phenotype, _qcovar, _covar, _extract, thin_common ->
        thin: thin_common
        supplied: true
    }
    LDAK_THINCOMMON(ch_fit_branches.thin.map { _fit_id, meta, bed, bim, fam, _phenotype, _qcovar, _covar, _extract, _thin_common -> [meta, bed, bim, fam] })

    def ch_thinned_fits = ch_fit_branches.thin
        .map { fit_id, meta, bed, bim, fam, phenotype, qcovar, covar, _extract, _thin_common -> [meta, fit_id, bed, bim, fam, phenotype, qcovar, covar] }
        .join(LDAK_THINCOMMON.out.predictors, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, fit_id, bed, bim, fam, phenotype, qcovar, covar, extract -> [fit_id, meta, bed, bim, fam, phenotype, qcovar, covar, extract] }
    def ch_fit_inputs = ch_fit_branches.supplied
        .map { fit_id, meta, bed, bim, fam, phenotype, qcovar, covar, extract, _thin_common -> [fit_id, meta, bed, bim, fam, phenotype, qcovar, covar, extract] }
        .mix(ch_thinned_fits)
    def ch_step1 = ch_fit_inputs.multiMap { _fit_id, meta, bed, bim, fam, phenotype, qcovar, covar, extract ->
        genotypes: [meta, bed, bim, fam]
        phenotype: [meta, phenotype]
        qcovar: [meta, qcovar]
        covar: [meta, covar]
        extract: [meta, extract]
    }
    LDAK_KVIKSTEP1(ch_step1.genotypes, ch_step1.phenotype, ch_step1.qcovar, ch_step1.covar, ch_step1.extract)
    def ch_fitted_predictions = LDAK_KVIKSTEP1.out.predictions
        .map { meta, root, details, prs -> [meta.id, root, details, prs] }
        .join(ch_fit_metadata, failOnDuplicate: true, failOnMismatch: true)
        .map { fit_id, root, details, prs, meta -> [fit_id, meta, root, details, prs] }
    def ch_predictors = LDAK_THINCOMMON.out.predictors
        .map { meta, extract -> [meta.id, extract] }
        .combine(ch_fit_metadata, by: 0)
        .map { fit_id, extract, meta -> [fit_id, meta, extract] }
    def ch_progress = LDAK_THINCOMMON.out.progress
        .map { meta, progress -> [meta.id, progress] }
        .combine(ch_fit_metadata, by: 0)
        .map { fit_id, progress, meta -> [fit_id, meta, progress] }
    def ch_effects = LDAK_KVIKSTEP1.out.effects
        .map { meta, effects -> [meta.id, effects] }
        .combine(ch_fit_metadata, by: 0)
        .map { fit_id, effects, meta -> [fit_id, meta, effects] }
    def ch_fit_log = LDAK_KVIKSTEP1.out.log
        .map { meta, log -> [meta.id, log] }
        .combine(ch_fit_metadata, by: 0)
        .map { fit_id, log, meta -> [fit_id, meta, log] }
    def ch_all_predictions = ch_fitted_predictions.mix(ch_predictions)
    def ch_step2 = ch_tests
        .combine(ch_all_predictions, by: 0)
        .multiMap { _fit_id, meta, bed, bim, fam, phenotype, qcovar, covar, keep, fit_meta, root, details, prs ->
            genotypes: [meta, bed, bim, fam]
            phenotype: [meta, phenotype]
            predictions: [fit_meta, root, details, prs]
            qcovar: [meta, qcovar]
            covar: [meta, covar]
            keep: [meta, keep]
        }
    LDAK_KVIKSTEP2(ch_step2.genotypes, ch_step2.phenotype, ch_step2.predictions, ch_step2.qcovar, ch_step2.covar, ch_step2.keep)

    emit:
    predictions = ch_all_predictions // channel: [ val(fit_id), val(meta), path(root), path(details), path(prs) ]
    predictors  = ch_predictors // channel: [ val(fit_id), val(meta), path(extract) ], generated predictor sets only
    progress    = ch_progress // channel: [ val(fit_id), val(meta), path(progress) ]
    results     = LDAK_KVIKSTEP2.out.results // channel: [ val(meta), path(assoc) ], native files per testing event
    summaries   = LDAK_KVIKSTEP2.out.summaries // channel: [ val(meta), path(summaries) ]
    pvalues     = LDAK_KVIKSTEP2.out.pvalues // channel: [ val(meta), path(pvalues) ]
    effects     = ch_effects // channel: [ val(fit_id), val(meta), path(effects) ]
    fit_log     = ch_fit_log // channel: [ val(fit_id), val(meta), path(log) ]
    test_log    = LDAK_KVIKSTEP2.out.log // channel: [ val(meta), path(log) ]
}
