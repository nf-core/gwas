// Route oriented relationship requests to shared sample/design preparation and MPH REML fitting.
// The pipeline supplies the shared pair tables, selects the matrix family, and preserves request identity.
// Trait orientation, missing-value spellings, native options, and publication remain caller-owned.

// SUBWORKFLOW: Pipeline-owned composition of local components
include { FIT_MPH_REML                    } from '../fit_mph_reml/main'

// FUNCTION: Local to the pipeline
include { getMethodTokensWithCapabilities } from '../validate_gwas_input/method_registry'

workflow ROUTE_MPH_BIVARIATE_RELATIONSHIPS {
    take:
    ch_relationships // channel: [ val(meta), [ path(genotype_file), ... ], path(pair_quant_covariates), path(pair_cat_covariates) ], one validated row per selected method per relationship
    ch_prepared_pairs // channel: [ val(relationship_meta), path(phenotype), path(quant_covariates), path(cat_covariates) ], headerless, one per relationship
    ch_named_covariates // channel: [ val(relationship_meta), path(quant_covariates), path(cat_covariates) ], headered, one per relationship; [] when absent
    ch_plink1_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam) ], relationship-scoped PLINK 1 bundles
    ch_mph_dense_matrices // channel: [ val(meta), val(matrix_identity), [ path(grm_file), ... ] ], relationship-scoped
    ch_mph_ldms_matrices // channel: [ val(meta), val(matrix_identity), [ path(grm_file), ... ], val(grm_prefixes) ], relationship-scoped

    main:

    // Selection is a registry query rather than a token literal: a request reaches this controller because its
    // pairwise entry declares the MPH option family, and it takes the dense or the stratified matrix stream
    // because of the matrix kind that entry declares.
    def mph_pair_tokens = getMethodTokensWithCapabilities([domain: 'pairwise', endpoint_domain: 'analysis', option_family: 'mph'])

    def ch_mph_requests = ch_relationships
        .filter { meta, _genotype_files, _pair_quant_covariates, _pair_cat_covariates -> meta.method in mph_pair_tokens }
        .map { meta, _genotype_files, _pair_quant_covariates, _pair_cat_covariates -> [meta.request_id, meta] }

    // One matrix record per request, whichever family produced it. The ordered prefix list is the family's
    // component order for the stratified case and the single bundle basename for the dense one, because
    // `--grm_list` names prefixes rather than paths and its order decides the result's row order.
    def ch_matrix_records = ch_mph_dense_matrices
        .map { meta, matrix_identity, grm_files -> [meta.request_id, matrix_identity, grm_files, [grmPrefix(grm_files)]] }
        .mix(
            ch_mph_ldms_matrices.map { meta, matrix_identity, grm_files, grm_prefixes -> [meta.request_id, matrix_identity, grm_files, grm_prefixes] }
        )

    // Both prepared streams carry exactly one record per relationship, so this is a guarded one-to-one join.
    // The fan-out to the several methods that may select one relationship is the `combine` below, which is
    // genuinely one-to-many and is the same seam the GCTA pair controller uses.
    def ch_prepared_relationships = ch_prepared_pairs
        .map { relationship_meta, phenotype, _quant_covariates, _cat_covariates -> [relationship_meta.relationship_id, phenotype] }
        .join(
            ch_named_covariates.map { relationship_meta, quant_covariates, cat_covariates ->
                [relationship_meta.relationship_id, quant_covariates ?: [], cat_covariates ?: []]
            },
            failOnDuplicate: true,
            failOnMismatch: true,
        )

    // Every other seam is a guarded join, so a request with no matrix, a matrix with no request and a request
    // with no PLINK 1 bundle each fail the run instead of vanishing.
    def ch_requests = ch_mph_requests
        .join(ch_matrix_records, failOnDuplicate: true, failOnMismatch: true)
        .map { request_id, meta, matrix_identity, grm_files, grm_prefixes ->
            [meta.relationship_id, request_id, resolveMphPairRouteMeta(meta, matrix_identity), grm_files, grm_prefixes]
        }
        .combine(ch_prepared_relationships, by: 0)
        .map { _relationship_id, request_id, route_meta, grm_files, grm_prefixes, phenotype, quant_covariates, cat_covariates ->
            [request_id, route_meta, grm_files, grm_prefixes, phenotype, quant_covariates, cat_covariates]
        }
        .join(
            ch_plink1_genotypes.filter { meta, _bed, _bim, _fam -> meta.method in mph_pair_tokens }.map { meta, _bed, _bim, fam -> [meta.request_id, fam] },
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { _request_id, route_meta, grm_files, grm_prefixes, phenotype, quant_covariates, cat_covariates, fam ->
            [route_meta, grm_files, grm_prefixes, phenotype, quant_covariates, cat_covariates, fam]
        }

    // The caller resolves matrix order and the accepted missing-value spelling before crossing the
    // portable preparation-and-fit boundary. Scientific method selection stays in this route.
    def ch_fit_inputs = ch_requests.map { route_meta, grm_files, grm_prefixes, phenotype, quant_covariates, cat_covariates, fam ->
        def ordered_grm_files = grm_prefixes.collectMany { grm_prefix ->
            grm_files.findAll { grm_file -> grm_file.name in ["${grm_prefix}.grm.bin".toString(), "${grm_prefix}.grm.iid".toString()] }
        }
        [route_meta, ordered_grm_files, phenotype, quant_covariates, cat_covariates, fam, mphTraitNames(route_meta), ['', 'na', 'nan', '-9']]
    }
    FIT_MPH_REML(ch_fit_inputs)

    emit:
    bivariate_results      = FIT_MPH_REML.out.variance_components.map { meta, variance_components -> [stripMphPairRouteState(meta), variance_components] } // channel: [ val(meta), path(mq.vc.csv) ], one per request
    bivariate_correlations = FIT_MPH_REML.out.correlations.map { meta, correlations -> [stripMphPairRouteState(meta), correlations] } // channel: [ val(meta), path(mq.cor.csv) ]
    bivariate_fixed        = FIT_MPH_REML.out.fixed_effects.map { meta, fixed_effects -> [stripMphPairRouteState(meta), fixed_effects] } // channel: [ val(meta), path(mq.blue.csv) ]
    bivariate_iterations   = FIT_MPH_REML.out.iterations.map { meta, iterations -> [stripMphPairRouteState(meta), iterations] } // channel: [ val(meta), path(mq.iter.csv) ]
    bivariate_log          = FIT_MPH_REML.out.native_log.map { meta, native_log -> [stripMphPairRouteState(meta), native_log] } // channel: [ val(meta), path(log) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Fold the matrix reuse identity and the MPH execution state into the request record. The request identity
// stays the primary metadata ID throughout execution: MPH addresses its matrix family through the `grm_list`
// it writes from the ordered prefixes, so nothing here needs the identity rewrite the GCTA LDMS route performs.
def resolveMphPairRouteMeta(pair_meta, matrix_identity) {
    return pair_meta + [
        matrix_key: matrix_identity.key,
        mph_estimator: pair_meta.method,
        mph_effective: buildMphEffectiveSettings(pair_meta.request_options.mph, pair_meta.native_args),
    ]
}

// The declared orientation of the pair, in manifest order and never sorted. This list becomes MPH's
// `--trait_names` and therefore fixes which endpoint is the left one in every published record.
def mphTraitNames(pair_meta) {
    return [pair_meta.left_trait_id, pair_meta.right_trait_id]
}

// The staged bundle's prefix, taken from the member whose name defines it. MPH's `--grm_list` names prefixes
// rather than paths, so this is the value the fit is addressed by.
def grmPrefix(grm_files) {
    def grm_iid_name = grm_files.find { grm_file -> grm_file.name.endsWith('.grm.iid') }.name
    return grm_iid_name.substring(0, grm_iid_name.length() - '.grm.iid'.length())
}

// Render the curated options and validated native additions for the configured MPH argument closure. Each
// fragment is rendered to a plain String: a GString is not equal to the String the task cache restores it as,
// so one left in metadata makes a resumed run's join by that metadata mismatch.
def buildMphEffectiveSettings(options, native_args) {
    return [native_arguments: [
        options.iterations != null ? "--num_iterations ${options.iterations}".toString() : null,
        options.tolerance != null ? "--tolerance ${options.tolerance}".toString() : null,
        options.random_vectors != null ? "--num_random_vectors ${options.random_vectors}".toString() : null,
        options.seed != null ? "--seed ${options.seed}".toString() : null,
        options.save_memory ? '--save_memory' : null,
    ].findAll { argument -> argument != null } + (native_args ?: []).collect { argument -> argument.toString() }]
}

// The route-local keys exist so the atoms, the adapters and the configured publication closures can see the
// selected token and its resolved settings. They are stripped before emission, so a consumer receives the
// request identity it supplied and the route cannot leak its own bookkeeping.
def stripMphPairRouteState(meta) {
    return meta.findAll { name, _value -> !(name in ['matrix_key', 'mph_estimator', 'mph_effective']) }
}
