// Route unary MPH requests to shared sample/design preparation and REML fitting.
// The pipeline owns estimator selection, matrix reuse, missing-value spellings, and result attribution.
// Native CSV preparation and fitting are composed by FIT_MPH_REML.

// SUBWORKFLOW: Pipeline-owned composition of local components
include { FIT_MPH_REML                    } from '../fit_mph_reml/main'

// FUNCTION: Local to the pipeline
include { getMethodTokensWithCapabilities } from '../validate_gwas_input/method_registry'

workflow ROUTE_MPH_HERITABILITY {
    take:
    ch_mph_dense_matrices // channel: [ val(meta), val(matrix_identity), [ path(grm_file), ... ] ]
    ch_mph_ldms_matrices // channel: [ val(meta), val(matrix_identity), [ path(grm_file), ... ], val(grm_prefixes) ]
    ch_plink1_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam) ], unary analyses needing PLINK 1
    ch_headerless_phenotypes // channel: [ val(meta), path(phenotype), path(quant_covariates), path(cat_covariates) ], covariates unused here
    ch_covariate_tables // channel: [ val(meta), path(quant_covariates), path(cat_covariates) ], headered; [] when absent

    main:

    // The estimator token comes from the registry query whose matrix kind produced the stream, so neither
    // token appears as a literal anywhere in this route.
    def ch_matrices = ch_mph_dense_matrices
        .map { meta, identity, grm_files -> [meta, identity, grm_files, [grmPrefix(grm_files)], resolveMphEstimator('mph_dense')] }
        .mix(
            ch_mph_ldms_matrices.map { meta, identity, grm_files, grm_prefixes -> [meta, identity, grm_files, grm_prefixes, resolveMphEstimator('mph_ldms')] }
        )

    // The matrix-by-phenotype seam is genuinely one-to-many -- one analysis may hold a dense record and an
    // LDMS record -- so it stays the `combine(by: 0)` the GRM heritability route uses. Every seam after it is
    // one-to-one on [analysis, estimator] and is a guarded join instead, because `combine` supports neither
    // `failOnMismatch` nor `failOnDuplicate` and would drop an unmatched analysis in silence.
    def ch_route_requests = ch_matrices
        .combine(ch_headerless_phenotypes.map { meta, phenotype, _quant_covariates, _cat_covariates -> [meta, phenotype] }, by: 0)
        .map { meta, _identity, grm_files, grm_prefixes, estimator, phenotype ->
            def route_meta = meta + [
                mph_estimator: estimator,
                mph_effective: buildMphEffectiveSettings(meta.method_options.mph),
            ]
            [routeKey(route_meta), route_meta, grm_files, grm_prefixes, phenotype]
        }

    // Both ancillary streams carry every analysis, LDAK and GCTA rows included, so each is first narrowed to
    // the rows that actually select an MPH token and fanned onto the same [analysis, estimator] key. Only then
    // is the join one-to-one in both directions and `failOnMismatch` meaningful: a selected MPH analysis with
    // no PLINK 1 bundle and a matrix record with no analysis both fail the run rather than vanishing.
    def ch_requests = ch_route_requests
        .join(
            ch_covariate_tables.flatMap { meta, quant_covariates, cat_covariates ->
                selectedMphEstimators(meta).collect { estimator -> [[meta.id, estimator], quant_covariates, cat_covariates] }
            },
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_plink1_genotypes.flatMap { meta, _bed, _bim, fam ->
                selectedMphEstimators(meta).collect { estimator -> [[meta.id, estimator], fam] }
            },
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { _key, route_meta, grm_files, grm_prefixes, phenotype, quant_covariates, cat_covariates, fam ->
            [route_meta, grm_files, grm_prefixes, phenotype, quant_covariates, cat_covariates, fam]
        }

    // The caller resolves matrix order, the fit's trait name and the accepted missing-value spelling before
    // crossing the portable preparation-and-fit boundary. A unary fit names the analysis's own declared trait,
    // so the list is built here rather than carried as a metadata key. Method selection stays in this route.
    def ch_fit_inputs = ch_requests.map { route_meta, grm_files, grm_prefixes, phenotype, quant_covariates, cat_covariates, fam ->
        def ordered_grm_files = grm_prefixes.collectMany { grm_prefix ->
            grm_files.findAll { grm_file -> grm_file.name in ["${grm_prefix}.grm.bin".toString(), "${grm_prefix}.grm.iid".toString()] }
        }
        [route_meta, ordered_grm_files, phenotype, quant_covariates, cat_covariates, fam, [route_meta.trait], ['', 'na', 'nan', '-9']]
    }
    FIT_MPH_REML(ch_fit_inputs)

    emit:
    mph_results    = FIT_MPH_REML.out.variance_components.map { meta, variance_components -> [stripMphRouteState(meta), variance_components] } // channel: [ val(meta), path(mq.vc.csv) ], one per analysis and selected MPH method
    mph_fixed      = FIT_MPH_REML.out.fixed_effects.map { meta, fixed_effects -> [stripMphRouteState(meta), fixed_effects] } // channel: [ val(meta), path(mq.blue.csv) ]
    mph_iterations = FIT_MPH_REML.out.iterations.map { meta, iterations -> [stripMphRouteState(meta), iterations] } // channel: [ val(meta), path(mq.iter.csv) ]
    mph_log        = FIT_MPH_REML.out.native_log.map { meta, native_log -> [stripMphRouteState(meta), native_log] } // channel: [ val(meta), path(log) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// The token whose matrix kind produced a stream. A registry query rather than a literal, and singular by
// assertion: two heritability tokens sharing one matrix kind would make the estimator of a result ambiguous.
def resolveMphEstimator(matrix_kind) {
    def tokens = getMethodTokensWithCapabilities([domain: 'heritability', option_family: 'mph', matrix_kind: matrix_kind])
    if (tokens.size() != 1) {
        error("[nf-core/gwas] ERROR: matrix kind '${matrix_kind}' resolves to ${tokens.size()} MPH heritability method(s) ${tokens}, expected exactly one")
    }
    return tokens.first()
}

// The row's own selected MPH tokens, used to fan the ancillary streams onto the [analysis, estimator] key the
// guarded joins are keyed by.
def selectedMphEstimators(meta) {
    def mph_methods = getMethodTokensWithCapabilities([domain: 'heritability', option_family: 'mph'])
    return meta.heritability_methods.findAll { method -> method in mph_methods }
}

def routeKey(meta) {
    return [meta.id, meta.mph_estimator]
}

// The staged bundle's prefix, taken from the member whose name defines it. MPH's `--grm_list` names prefixes
// rather than paths, so this is the value the fit is addressed by.
def grmPrefix(grm_files) {
    def grm_iid_name = grm_files.find { grm_file -> grm_file.name.endsWith('.grm.iid') }.name
    return grm_iid_name.substring(0, grm_iid_name.length() - '.grm.iid'.length())
}

// Render the curated options for the configured MPH argument closure. Each fragment is rendered to a plain
// String: a GString is not equal to the String the task cache restores it as, so one left in metadata makes a
// resumed run's join by that metadata mismatch.
def buildMphEffectiveSettings(options) {
    return [native_arguments: [
        options.iterations != null ? "--num_iterations ${options.iterations}".toString() : null,
        options.tolerance != null ? "--tolerance ${options.tolerance}".toString() : null,
        options.random_vectors != null ? "--num_random_vectors ${options.random_vectors}".toString() : null,
        options.seed != null ? "--seed ${options.seed}".toString() : null,
        options.save_memory ? '--save_memory' : null,
    ].findAll { argument -> argument != null }]
}

// The route-local keys exist so the atoms, the adapters and the configured publication closures can see the
// selected token and its resolved settings. Everything but the token is stripped before emission, so a
// consumer receives the focal analysis identity it supplied and the route cannot leak its own bookkeeping.
def stripMphRouteState(meta) {
    return meta.findAll { name, _value -> !(name in ['mph_effective']) }
}
