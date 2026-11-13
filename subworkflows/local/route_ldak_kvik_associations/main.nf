// Route nf-core/gwas analysis records through LDAK-KVIK with shared genotype-level predictor resources.
// This is pipeline-specific relational-input policy, not an nf-core/modules submission candidate.
// Every constituent process reports directly to the run-wide versions topic, so this subworkflow emits no versions.

// MODULES: Upstream-ready components used inside a pipeline-local route
include { LDAK_THINCOMMON                 } from '../../../modules/local/ldak/thincommon/main'
include { PREPARE_LDAK_KVIK_HARMONISATION } from '../../../modules/local/prepare_ldak_kvik_harmonisation/main'
include { PLINK_ASSOCIATION_LDAK_KVIK     } from '../plink_association_ldak_kvik/main'

// FUNCTION: Local to the pipeline
include { digestFileBytes                 } from '../utils_nfcore_gwas_pipeline'
include { buildScientificArtifactKey      } from '../utils_nfcore_gwas_pipeline'

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
    // Predictor artifacts are resolved before fitting. Thin-common requests collapse on genotype
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
            // order would then invalidate all of them.
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

    // The manifest declares one fit per analysis. Predictor resources remain shared by genotype view.
    def ch_requests = ch_predictor_requests
        .combine(ch_resolved_predictors, by: 0)
        .map { _predictor_artifact_key, meta, _view_key, bed, bim, fam, _predictor_artifact, extract ->
            [meta.id, meta, bed, bim, fam, extract]
        }
        .join(ch_phenotypes.map { meta, phenotype, qcovar, covar -> [meta.id, phenotype, qcovar, covar] }, failOnDuplicate: true, failOnMismatch: true)
    def ch_fits = ch_requests.map { fit_id, meta, bed, bim, fam, extract, phenotype, qcovar, covar ->
        // Fit metadata excludes Step 2-only policy so that changing a testing sample subset preserves the fit.
        [fit_id, [id: meta.id, is_binary: meta.is_binary], bed, bim, fam, phenotype, qcovar, covar, extract, false]
    }
    def ch_tests = ch_requests.map { fit_id, meta, bed, bim, fam, _extract, phenotype, qcovar, covar ->
        [fit_id, meta, bed, bim, fam, phenotype, qcovar, covar, meta.method_options.ldak.kvik_step2_keep]
    }
    PLINK_ASSOCIATION_LDAK_KVIK(ch_fits, channel.empty(), ch_tests)

    // The pipeline's one test bundle per analysis makes this an exact join; public P8 may emit many shards.
    def ch_harmonisation = PLINK_ASSOCIATION_LDAK_KVIK.out.results.join(PLINK_ASSOCIATION_LDAK_KVIK.out.summaries, failOnDuplicate: true, failOnMismatch: true)
    PREPARE_LDAK_KVIK_HARMONISATION(ch_harmonisation)

    emit:
    harmonisation_input = PREPARE_LDAK_KVIK_HARMONISATION.out.harmonisation_input // channel: [ val(meta), path(tsv) ]
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
