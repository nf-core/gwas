// Route nf-core/gwas analysis records through LDAK-KVIK while reusing scientifically identical Step 1 fits.
// This is pipeline-specific relational-input and publication policy, not an nf-core/modules submission candidate.
// Every constituent process reports directly to the run-wide versions topic, so this subworkflow emits no versions.

// MODULES: Upstream-ready components used inside a pipeline-local route
include { LDAK_THINCOMMON                 } from '../../../modules/local/ldak/thincommon/main'
include { LDAK_KVIKSTEP1                  } from '../../../modules/local/ldak/kvikstep1/main'
include { LDAK_KVIKSTEP2                  } from '../../../modules/local/ldak/kvikstep2/main'
include { ATTRIBUTE_LDAK_KVIK_PREDICTIONS } from '../../../modules/local/attribute_ldak_kvik_predictions/main'

// FUNCTION: Local to the pipeline
include { digestFileBytes                 } from '../utils_nfcore_gwas_pipeline'
include { buildCanonicalPredictionKey     } from '../utils_nfcore_gwas_pipeline'

workflow ROUTE_LDAK_KVIK_ASSOCIATIONS {
    take:
    ch_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam) ], once per analysis
    ch_phenotypes // channel: [ val(meta2), path(phenotype), path(quant_covariates), path(cat_covariates) ], optional files are []
    ch_extract_policy // channel: [ val(meta3), path(extract), val(subset_policy) ], extract is [] for all/thin_common

    main:
    // The reuse key is a private routing value rather than a custom meta field. One canonical request
    // drives each fit; the original analysis metadata stays beside every consumer and is restored below.
    def ch_requests = ch_genotypes
        .map { meta, bed, bim, fam -> [meta.id, meta, bed, bim, fam] }
        .join(ch_phenotypes.map { meta, phenotype, quant_covariates, cat_covariates -> [meta.id, phenotype, quant_covariates, cat_covariates] }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_extract_policy.map { meta, extract, subset_policy -> [meta.id, extract, subset_policy] }, failOnDuplicate: true, failOnMismatch: true)
        .map { _analysis_id, meta, bed, bim, fam, phenotype, quant_covariates, cat_covariates, extract, subset_policy ->
            def prediction_key = buildKvikPredictionKey(meta, phenotype, quant_covariates, cat_covariates, subset_policy, extract)
            def fit_meta = meta + [id: "${meta.cohort}.ldak_kvik.${prediction_key}"]
            [prediction_key, meta, fit_meta, bed, bim, fam, phenotype, quant_covariates, cat_covariates, extract, subset_policy]
        }

    def ch_fit_requests = ch_requests.unique { prediction_key, _meta, _fit_meta, _bed, _bim, _fam, _phenotype, _quant_covariates, _cat_covariates, _extract, _subset_policy -> prediction_key }

    def ch_thin_common_inputs = ch_fit_requests.filter { _prediction_key, _meta, _fit_meta, _bed, _bim, _fam, _phenotype, _quant_covariates, _cat_covariates, _extract, subset_policy -> subset_policy == 'thin_common' }
    LDAK_THINCOMMON(ch_thin_common_inputs.map { _prediction_key, _meta, fit_meta, bed, bim, fam, _phenotype, _quant_covariates, _cat_covariates, _extract, _subset_policy -> [fit_meta, bed, bim, fam] })

    def ch_direct_step1_inputs = ch_fit_requests
        .filter { _prediction_key, _meta, _fit_meta, _bed, _bim, _fam, _phenotype, _quant_covariates, _cat_covariates, _extract, subset_policy -> subset_policy in ['all', 'provided'] }
        .map { _prediction_key, _meta, fit_meta, bed, bim, fam, phenotype, quant_covariates, cat_covariates, extract, _subset_policy ->
            [[fit_meta, bed, bim, fam], [fit_meta, phenotype, fit_meta.is_binary], [fit_meta, quant_covariates], [fit_meta, cat_covariates], [fit_meta, extract]]
        }

    def ch_thin_step1_inputs = ch_thin_common_inputs
        .map { _prediction_key, _meta, fit_meta, bed, bim, fam, phenotype, quant_covariates, cat_covariates, _extract, _subset_policy ->
            [fit_meta.id, [fit_meta, bed, bim, fam], [fit_meta, phenotype, fit_meta.is_binary], [fit_meta, quant_covariates], [fit_meta, cat_covariates]]
        }
        .join(LDAK_THINCOMMON.out.predictors.map { meta, extract -> [meta.id, [meta, extract]] }, failOnDuplicate: true, failOnMismatch: true)
        .map { _fit_id, genotypes, pheno, qcovar, covar, extract -> [genotypes, pheno, qcovar, covar, extract] }

    def ch_step1 = ch_direct_step1_inputs
        .mix(ch_thin_step1_inputs)
        .multiMap { genotypes, pheno, qcovar, covar, extract ->
            genotypes: genotypes
            pheno: pheno
            qcovar: qcovar
            covar: covar
            extract: extract
        }

    LDAK_KVIKSTEP1(ch_step1.genotypes, ch_step1.pheno, ch_step1.qcovar, ch_step1.covar, ch_step1.extract)

    // Reattach the private prediction key through the synthetic fit id, then fan the fitted bundle
    // out to every consuming analysis while restoring its original focal metadata for Step 2.
    def ch_fit_keys = ch_fit_requests.map { prediction_key, _meta, fit_meta, _bed, _bim, _fam, _phenotype, _quant_covariates, _cat_covariates, _extract, _subset_policy -> [fit_meta.id, prediction_key] }
    def ch_prediction_bundles = LDAK_KVIKSTEP1.out.predictions
        .map { fit_meta, root, loco_details, loco_prs -> [fit_meta.id, fit_meta, root, loco_details, loco_prs] }
        .join(ch_fit_keys, failOnDuplicate: true, failOnMismatch: true)
        .map { _fit_id, fit_meta, root, loco_details, loco_prs, prediction_key -> [prediction_key, [fit_meta, root, loco_details, loco_prs]] }

    def ch_attributed_predictions = ch_requests
        .map { prediction_key, meta, _fit_meta, _bed, _bim, _fam, _phenotype, _quant_covariates, _cat_covariates, _extract, _subset_policy -> [prediction_key, meta.id, meta] }
        .unique { _prediction_key, analysis_id, _meta -> analysis_id }
        .combine(ch_prediction_bundles, by: 0)
        .map { _prediction_key, _analysis_id, meta, predictions -> [meta, predictions[1], predictions[2], predictions[3]] }

    ATTRIBUTE_LDAK_KVIK_PREDICTIONS(ch_attributed_predictions)

    def ch_step2 = ch_requests
        .combine(ch_prediction_bundles, by: 0)
        .multiMap { _prediction_key, meta, _fit_meta, bed, bim, fam, phenotype, quant_covariates, cat_covariates, _extract, _subset_policy, predictions ->
            genotypes: [meta, bed, bim, fam]
            pheno: [meta, phenotype]
            predictions: predictions
            qcovar: [meta, quant_covariates]
            covar: [meta, cat_covariates]
            keep: [meta, []]
        }

    LDAK_KVIKSTEP2(ch_step2.genotypes, ch_step2.pheno, ch_step2.predictions, ch_step2.qcovar, ch_step2.covar, ch_step2.keep)

    emit:
    results             = LDAK_KVIKSTEP2.out.results // channel: [ val(meta), path(assoc) ]
    harmonisation_input = LDAK_KVIKSTEP2.out.harmonisation_input // channel: [ val(meta), path(tsv) ]
    summaries           = LDAK_KVIKSTEP2.out.summaries // channel: [ val(meta), path(summaries) ]
    pvalues             = LDAK_KVIKSTEP2.out.pvalues // channel: [ val(meta), path(pvalues) ], optional
    predictions         = LDAK_KVIKSTEP1.out.predictions // channel: [ val(meta), path(root), path(loco_details), path(loco_prs) ], once per shared fit
    effects             = LDAK_KVIKSTEP1.out.effects // channel: [ val(meta), path(effects) ], optional per shared fit
    progress            = LDAK_THINCOMMON.out.progress // channel: [ val(meta), path(progress) ], only for thin_common fits
    logs                = LDAK_KVIKSTEP1.out.log.mix(LDAK_KVIKSTEP2.out.log) // channel: [ val(meta), path(log) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// LDAK-KVIK Step 1 reuse additionally depends on predictor policy and optional predictor-list bytes.
def buildKvikPredictionKey(meta, phenotype, quant_covariates, cat_covariates, subset_policy, predictor_extract) {
    def identity = [
        cohort: meta.cohort,
        trait: meta.trait,
        is_binary: meta.is_binary,
        phenotype: getPredictionInputIdentity(phenotype),
        quant_covariates: getPredictionInputIdentity(quant_covariates),
        cat_covariates: getPredictionInputIdentity(cat_covariates),
        subset_policy: subset_policy,
        predictor_extract: getPredictionInputIdentity(predictor_extract),
    ]
    return buildCanonicalPredictionKey(identity)
}

def getPredictionInputIdentity(input_file) {
    return input_file ? digestFileBytes(input_file) : 'absent'
}
