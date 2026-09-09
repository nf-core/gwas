// Route nf-core/gwas analysis records through LDAK-KVIK while reusing scientifically identical Step 1 fits.
// This is pipeline-specific relational-input policy, not an nf-core/modules submission candidate.
// Every constituent process reports directly to the run-wide versions topic, so this subworkflow emits no versions.

// MODULES: Upstream-ready components used inside a pipeline-local route
include { LDAK_THINCOMMON             } from '../../../modules/local/ldak/thincommon/main'
include { LDAK_KVIKSTEP1              } from '../../../modules/local/ldak/kvikstep1/main'
include { LDAK_KVIKSTEP2              } from '../../../modules/local/ldak/kvikstep2/main'

// FUNCTION: Local to the pipeline
include { digestFileBytes             } from '../utils_nfcore_gwas_pipeline'
include { buildCanonicalPredictionKey } from '../utils_nfcore_gwas_pipeline'
include { buildScientificArtifactKey  } from '../utils_nfcore_gwas_pipeline'

workflow ROUTE_LDAK_KVIK_ASSOCIATIONS {
    take:
    ch_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam), val(view_key) ], once per analysis
    ch_phenotypes // channel: [ val(meta2), path(phenotype), path(quant_covariates), path(cat_covariates) ], optional files are []
    ch_extract_policy // channel: [ val(meta3), path(extract), val(subset_policy) ], extract is [] for all/thin_common

    main:
    // A supplied predictor resource is hashed once per distinct declared path before it reaches the
    // phenotype-specific request graph. Byte-identical files still share the resulting artifact identity.
    def ch_provided_predictor_identities = ch_extract_policy
        .filter { _meta, _extract, subset_policy -> subset_policy == 'provided' }
        .map { _meta, extract, _subset_policy -> [extract.toString(), extract] }
        .unique { resource_path, _extract -> resource_path }
        .map { resource_path, extract -> [resource_path, buildProvidedPredictorArtifactKey(extract)] }

    def ch_predictor_policies = ch_extract_policy
        .filter { _meta, _extract, subset_policy -> subset_policy != 'provided' }
        .map { meta, extract, subset_policy -> [meta.id, extract, subset_policy, []] }

    def ch_provided_predictor_policies = ch_extract_policy
        .filter { _meta, _extract, subset_policy -> subset_policy == 'provided' }
        .map { meta, extract, subset_policy -> [extract.toString(), meta.id, extract, subset_policy] }
        .combine(ch_provided_predictor_identities, by: 0)
        .map { _resource_path, analysis_id, extract, subset_policy, predictor_artifact_key ->
            [analysis_id, extract, subset_policy, predictor_artifact_key]
        }
    ch_predictor_policies = ch_predictor_policies.mix(ch_provided_predictor_policies)

    // The PLINK 1 view key arrives from preparation as a tuple member. It is the cohort's immutable PLINK 1
    // view identity — the supplied bundle's own identity for a native PLINK 1 cohort, the hard-call
    // projection's for a PLINK 2 or VCF one — so this route infers nothing from a staged basename and
    // absorbs no backend policy.
    //
    // Predictor artifacts are resolved before Step 1 identity. Thin-common requests collapse on genotype
    // view plus the fixed native thinning contract; phenotype, covariates and focal analysis identity are not
    // present at this level.
    def ch_predictor_requests = ch_genotypes
        .map { meta, bed, bim, fam, view_key -> [meta.id, meta, view_key, bed, bim, fam] }
        .join(ch_predictor_policies, failOnDuplicate: true, failOnMismatch: true)
        .map { _analysis_id, meta, view_key, bed, bim, fam, extract, subset_policy, provided_artifact_key ->
            def predictor_artifact = buildKvikPredictorArtifact(view_key, subset_policy, extract, provided_artifact_key)
            [predictor_artifact.key, meta, view_key, bed, bim, fam, predictor_artifact]
        }

    def ch_thin_common_inputs = ch_predictor_requests
        .filter { _predictor_artifact_key, _meta, _view_key, _bed, _bim, _fam, predictor_artifact -> predictor_artifact.policy == 'thin_common' }
        .unique { predictor_artifact_key, _meta, _view_key, _bed, _bim, _fam, _predictor_artifact -> predictor_artifact_key }
        .map { predictor_artifact_key, _meta, _view_key, bed, bim, fam, _predictor_artifact ->
            // The predictor artifact key already folds the genotype view and the fixed thinning contract, so
            // it is unique on its own — and it must be the whole identity. The `unique` above collapses every
            // cohort sharing one PLINK 1 view onto one element, so naming the artifact after the surviving
            // element's cohort would put an arbitrary cohort into this task's hash, its output filenames and,
            // through the predictor file staged into Step 1, every fit that consumes it. A change in arrival
            // order would then invalidate all of them. This is the same reason `base.cohort` was removed from
            // the relatedness artifact identity, and the Step 1 fit meta below is cohort-free for it too.
            def thin_meta = [id: "ldak_thin_common.${predictor_artifact_key}", predictor_artifact_key: predictor_artifact_key]
            [thin_meta, bed, bim, fam]
        }
    LDAK_THINCOMMON(ch_thin_common_inputs)

    def ch_resolved_predictors = ch_predictor_requests
        .filter { _predictor_artifact_key, _meta, _view_key, _bed, _bim, _fam, predictor_artifact -> predictor_artifact.policy != 'thin_common' }
        .unique { predictor_artifact_key, _meta, _view_key, _bed, _bim, _fam, _predictor_artifact -> predictor_artifact_key }
        .map { predictor_artifact_key, _meta, _view_key, _bed, _bim, _fam, predictor_artifact ->
            [predictor_artifact_key, predictor_artifact.extract]
        }
    ch_resolved_predictors = ch_resolved_predictors.mix(
        LDAK_THINCOMMON.out.predictors.map { meta, extract -> [meta.predictor_artifact_key, extract] }
    )

    // Step 1 identity is computed only after the predictor resource has been resolved. It names the selected
    // artifact directly, so all, provided and thin-common policies cannot collapse onto one fit identity.
    def ch_resolved_predictor_requests = ch_predictor_requests
        .combine(ch_resolved_predictors, by: 0)
        .map { _predictor_artifact_key, meta, view_key, bed, bim, fam, predictor_artifact, resolved_extract ->
            [meta.id, meta, view_key, bed, bim, fam, predictor_artifact, resolved_extract]
        }

    def ch_resolved_requests = ch_resolved_predictor_requests
        .join(ch_phenotypes.map { meta, phenotype, quant_covariates, cat_covariates -> [meta.id, phenotype, quant_covariates, cat_covariates] }, failOnDuplicate: true, failOnMismatch: true)
        .map { _analysis_id, meta, view_key, bed, bim, fam, predictor_artifact, resolved_extract, phenotype, quant_covariates, cat_covariates ->
            def prediction_key = buildKvikPredictionKey(meta, view_key, phenotype, quant_covariates, cat_covariates, predictor_artifact.key)
            def fit_meta = [id: "ldak_kvik.${prediction_key}", is_binary: meta.is_binary]
            [prediction_key, meta, fit_meta, bed, bim, fam, phenotype, quant_covariates, cat_covariates, resolved_extract]
        }

    def ch_fit_requests = ch_resolved_requests
        .map { prediction_key, meta, fit_meta, bed, bim, fam, phenotype, quant_covariates, cat_covariates, extract ->
            [prediction_key, [prediction_key, meta, fit_meta, bed, bim, fam, phenotype, quant_covariates, cat_covariates, extract]]
        }
        .groupTuple()
        .map { _prediction_key, requests -> requests.sort { left, right -> left[1].id <=> right[1].id }.first() }

    def ch_step1 = ch_fit_requests.multiMap { _prediction_key, _meta, fit_meta, bed, bim, fam, phenotype, quant_covariates, cat_covariates, extract ->
        genotypes: [fit_meta, bed, bim, fam]
        pheno: [fit_meta, phenotype, fit_meta.is_binary]
        qcovar: [fit_meta, quant_covariates]
        covar: [fit_meta, cat_covariates]
        extract: [fit_meta, extract]
    }

    LDAK_KVIKSTEP1(ch_step1.genotypes, ch_step1.pheno, ch_step1.qcovar, ch_step1.covar, ch_step1.extract)

    // Reattach the private prediction key through the synthetic fit id, then fan the fitted bundle
    // out to every consuming analysis while restoring its original focal metadata for Step 2.
    def ch_fit_keys = ch_fit_requests.map { prediction_key, _meta, fit_meta, _bed, _bim, _fam, _phenotype, _quant_covariates, _cat_covariates, _extract -> [fit_meta.id, prediction_key] }
    def ch_prediction_bundles = LDAK_KVIKSTEP1.out.predictions
        .map { fit_meta, root, loco_details, loco_prs -> [fit_meta.id, fit_meta, root, loco_details, loco_prs] }
        .join(ch_fit_keys, failOnDuplicate: true, failOnMismatch: true)
        .map { _fit_id, fit_meta, root, loco_details, loco_prs, prediction_key -> [prediction_key, [fit_meta, root, loco_details, loco_prs]] }

    def ch_step2 = ch_resolved_requests
        .combine(ch_prediction_bundles, by: 0)
        .multiMap { _prediction_key, meta, _fit_meta, bed, bim, fam, phenotype, quant_covariates, cat_covariates, _extract, predictions ->
            genotypes: [meta, bed, bim, fam]
            pheno: [meta, phenotype]
            predictions: predictions
            qcovar: [meta, quant_covariates]
            covar: [meta, cat_covariates]
            keep: [meta, meta.method_options.ldak.kvik_step2_keep]
        }

    LDAK_KVIKSTEP2(ch_step2.genotypes, ch_step2.pheno, ch_step2.predictions, ch_step2.qcovar, ch_step2.covar, ch_step2.keep)

    emit:
    results             = LDAK_KVIKSTEP2.out.results // channel: [ val(meta), path(assoc) ]
    harmonisation_input = LDAK_KVIKSTEP2.out.harmonisation_input // channel: [ val(meta), path(tsv) ]
    progress            = LDAK_THINCOMMON.out.progress // channel: [ val(meta), path(progress) ], once per thin-common predictor artifact
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def buildThinCommonKey(view_key, effective_native_options = []) {
    return "thin_common.${buildScientificArtifactKey(
        [type: 'ldak_thin_common_predictors', view_key: view_key, adapter_contract: 'ldak6_6.1_thincommon_v1'],
        effective_native_options,
    )}"
}

def buildProvidedPredictorArtifactKey(predictor_extract) {
    return "provided.${digestFileBytes(predictor_extract)}"
}

def buildKvikPredictorArtifact(view_key, subset_policy, predictor_extract, provided_artifact_key) {
    if (subset_policy == 'all') {
        return [key: 'all_predictors', policy: subset_policy, extract: []]
    }
    if (subset_policy == 'provided') {
        return [key: provided_artifact_key, policy: subset_policy, extract: predictor_extract]
    }
    return [key: buildThinCommonKey(view_key), policy: subset_policy, extract: []]
}

// LDAK-KVIK Step 1 reuse depends on the resolved predictor artifact, not ownership of its construction.
def buildKvikPredictionKey(meta, view_key, phenotype, quant_covariates, cat_covariates, predictor_artifact_key) {
    def identity = [
        genotype_view: view_key,
        is_binary: meta.is_binary,
        phenotype: getPredictionInputIdentity(phenotype),
        quant_covariates: getPredictionInputIdentity(quant_covariates),
        cat_covariates: getPredictionInputIdentity(cat_covariates),
        predictor_artifact: predictor_artifact_key,
        adapter_contract: 'ldak6_6.1_kvikstep1_v1',
        native_options: '--mpheno 1',
    ]
    return buildCanonicalPredictionKey(identity)
}

def getPredictionInputIdentity(input_file) {
    return input_file ? digestFileBytes(input_file) : 'absent'
}
