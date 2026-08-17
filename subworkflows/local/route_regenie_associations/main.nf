// Route nf-core/gwas analysis records through REGENIE while reusing scientifically identical Step 1 fits.
// This is pipeline-specific relational-input and publication policy, not an nf-core/modules submission candidate.
// Every constituent process reports directly to the run-wide versions topic, so this subworkflow emits no versions.

// SUBWORKFLOW: Upstream-ready REGENIE composition used inside a pipeline-local route
include { PLINK_FIT_REGENIE             } from '../plink_fit_regenie/main'

// MODULE: Local to the pipeline
include { ATTRIBUTE_REGENIE_PREDICTIONS } from '../../../modules/local/attribute_regenie_predictions/main'

// MODULE: Installed directly from nf-core/modules
include { REGENIE_STEP2                 } from '../../../modules/nf-core/regenie/step2/main'

// FUNCTION: Local to the pipeline
include { digestFileBytes               } from '../utils_nfcore_gwas_pipeline'
include { buildCanonicalPredictionKey   } from '../utils_nfcore_gwas_pipeline'

workflow ROUTE_REGENIE_ASSOCIATIONS {
    take:
    ch_analyses // channel: [ val(meta), path(pgen), path(psam), path(pvar), path(phenotype), path(covariates) ], covariates is [] when absent
    step2_bsize // channel: val(step2_bsize)
    step1_mode // channel: val(step1_mode), 'standard' or 'chunked'
    step1_jobs // channel: val(step1_jobs), use null for standard mode

    main:

    if (step1_mode == 'chunked' && step1_jobs == null) {
        error("[nf-core/gwas] ERROR: --regenie_step1_jobs is required when --regenie_step1_mode is 'chunked'")
    }

    // The reuse key is a private routing value rather than a custom meta field. One canonical request
    // drives each fit; the original analysis metadata stays beside every consumer and is restored below.
    def ch_requests = ch_analyses.map { meta, pgen, psam, pvar, phenotype, covariates ->
        def step1_bsize = meta.method_options.regenie.step1_bsize
        def prediction_key = buildRegeniePredictionKey(meta, phenotype, covariates ?: [], step1_bsize)
        def fit_meta = meta + [id: "${meta.cohort}.regenie.${prediction_key}"]
        [prediction_key, meta, fit_meta, pgen, pvar, psam, phenotype, covariates ?: [], step1_bsize]
    }

    def ch_fit_requests = ch_requests.unique { prediction_key, _meta, _fit_meta, _pgen, _pvar, _psam, _phenotype, _covariates, _step1_bsize -> prediction_key }

    def ch_fit = ch_fit_requests.multiMap { _prediction_key, _meta, fit_meta, pgen, pvar, psam, phenotype, covariates, step1_bsize ->
        genotypes: [fit_meta, pgen, pvar, psam]
        pheno: [fit_meta, phenotype]
        covar: [fit_meta, covariates]
        bsize: [fit_meta, step1_bsize]
        mode: [fit_meta, step1_mode]
        jobs: [fit_meta, step1_mode == 'chunked' ? step1_jobs : []]
    }

    PLINK_FIT_REGENIE(ch_fit.genotypes, ch_fit.pheno, ch_fit.covar, ch_fit.bsize, ch_fit.mode, ch_fit.jobs)

    // PLINK_FIT_REGENIE returns only opaque fit metadata. Reattach the private key through the fit id
    // side channel, then combine one fitted bundle with every analysis request that consumes it.
    def ch_fit_keys = ch_fit_requests.map { prediction_key, _meta, fit_meta, _pgen, _pvar, _psam, _phenotype, _covariates, _step1_bsize ->
        [fit_meta.id, prediction_key]
    }

    def ch_prediction_bundles = PLINK_FIT_REGENIE.out.predictions
        .map { fit_meta, predictions -> [fit_meta.id, predictions] }
        .join(
            PLINK_FIT_REGENIE.out.loco.map { fit_meta, loco -> [fit_meta.id, loco] },
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(ch_fit_keys, failOnDuplicate: true, failOnMismatch: true)
        .map { _fit_id, predictions, loco, prediction_key -> [prediction_key, predictions, loco] }

    def ch_attributed_predictions = ch_requests
        .map { prediction_key, meta, _fit_meta, _pgen, _pvar, _psam, _phenotype, _covariates, _step1_bsize -> [prediction_key, meta.id, meta] }
        .unique { _prediction_key, analysis_id, _meta -> analysis_id }
        .combine(ch_prediction_bundles, by: 0)
        .map { _prediction_key, _analysis_id, meta, predictions, loco -> [meta, predictions, loco] }

    ATTRIBUTE_REGENIE_PREDICTIONS(ch_attributed_predictions)

    def ch_step2 = ch_requests
        .combine(ch_prediction_bundles, by: 0)
        .multiMap { _prediction_key, meta, _fit_meta, pgen, pvar, psam, phenotype, covariates, _step1_bsize, predictions, loco ->
            genotypes: [meta, pgen, pvar, psam]
            predictions: [meta, predictions, loco]
            pheno: [meta, phenotype]
            covar: [meta, covariates]
            bsize: step2_bsize
        }

    REGENIE_STEP2(ch_step2.genotypes, ch_step2.predictions, ch_step2.pheno, ch_step2.covar, ch_step2.bsize)

    emit:
    results     = REGENIE_STEP2.out.results // channel: [ val(meta), path(regenie_results) ]
    logs        = PLINK_FIT_REGENIE.out.logs.mix(REGENIE_STEP2.out.log) // channel: [ val(meta), path(log) ]
    predictions = ch_attributed_predictions.map { meta, predictions, _loco -> [meta, predictions] } // channel: [ val(meta), path(predictions) ]
    loco        = ch_attributed_predictions.map { meta, _predictions, loco -> [meta, loco] } // channel: [ val(meta), path(loco) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// REGENIE Step 1 reuse requires every cohort-defining and scientific input to agree. Execution-only
// controls are absent; Step 1 block size remains because it changes the fitted model.
def buildRegeniePredictionKey(meta, phenotype, covariates, step1_bsize) {
    def identity = [
        cohort: meta.cohort,
        trait: meta.trait,
        is_binary: meta.is_binary,
        phenotype: getPredictionInputIdentity(phenotype),
        covariates: getPredictionInputIdentity(covariates),
        step1_bsize: step1_bsize,
    ]
    return buildCanonicalPredictionKey(identity)
}

def getPredictionInputIdentity(input_file) {
    return input_file ? digestFileBytes(input_file) : 'absent'
}
