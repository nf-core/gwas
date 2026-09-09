// Route nf-core/gwas individual-level relationship requests through MPH bivariate REML on native MPH
// relatedness matrices, one-component and LD-by-MAF-stratified. A relationship is a declared, oriented pair of
// analysis units over one cohort. The ordered two-trait table and the relationship-owned covariates are
// prepared once on the spine and shared with the GCTA bivariate controller, so selecting MPH beside GCTA costs
// no second preparation and defines no second endpoint sample set.
//
// MPH matrices are a different on-disk format from GCTA's and are never interchanged: this controller consumes
// only the mph_dense and mph_ldms families and takes each request's stream from the registry's matrix kind,
// never from a method-name literal. It owns estimator selection, the IID-keyed serialisation of the shared
// pair and the effective settings the native call renders.
//
// Declared left/right orientation is written into `--trait_names` in that order. Native output labels are
// retained unchanged. This controller owns pipeline routing, reuse identity and publication policy.
// Every constituent process reports directly to the run-wide versions topic, so this subworkflow emits no
// versions, and it reads no params, no workflow and no projectDir.

// MODULE: Local to the pipeline
include { PREPARE_MPH_INPUTS              } from '../../../modules/local/prepare_mph_inputs/main'
include { MPH_REML                        } from '../../../modules/local/mph/reml/main'

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

    //
    // MODULE: Serialise the shared pair into the IID-keyed CSVs MPH consumes, in the matrix family's own order
    //
    def ch_serializer_inputs = ch_requests.multiMap { route_meta, grm_files, _grm_prefixes, phenotype, quant_covariates, cat_covariates, fam ->
        tables: [route_meta, phenotype, quant_covariates, cat_covariates]
        // Every component shares one matrix sample order; the FAM supplies the matching FID for each IID.
        identity: [route_meta, grm_files.find { grm_file -> grm_file.name.endsWith('.grm.iid') }, fam]
    }
    PREPARE_MPH_INPUTS(ch_serializer_inputs.tables, ch_serializer_inputs.identity)

    //
    // MODULE: Fit the two traits jointly over the matrix family
    //
    def ch_reml_inputs = ch_requests
        .map { route_meta, grm_files, grm_prefixes, _phenotype, _quant_covariates, _cat_covariates, _fam -> [route_meta.request_id, route_meta, grm_files, grm_prefixes] }
        .join(PREPARE_MPH_INPUTS.out.phenotype.map { meta, phenotype_csv -> [meta.request_id, phenotype_csv] }, failOnDuplicate: true, failOnMismatch: true)
        .join(PREPARE_MPH_INPUTS.out.covariates.map { meta, covariate_csv -> [meta.request_id, covariate_csv] }, failOnDuplicate: true, remainder: true)
        .filter { record -> record[1] != null }
        .multiMap { _request_id, route_meta, grm_files, grm_prefixes, phenotype_csv, covariate_csv ->
            grm: [route_meta, grm_files, grm_prefixes]
            // The declared left/right order fixes the native trait orientation.
            pheno: [route_meta, phenotype_csv, route_meta.mph_trait_names]
            // Read the encoded column names from the completed CSV.
            covar: [route_meta, covariate_csv ?: [], covariate_csv ? readCsvHeader(covariate_csv).drop(1) : []]
        }
    MPH_REML(ch_reml_inputs.grm, ch_reml_inputs.pheno, ch_reml_inputs.covar)

    emit:
    bivariate_results      = MPH_REML.out.variance_components.map { meta, variance_components -> [stripMphPairRouteState(meta), variance_components] } // channel: [ val(meta), path(mq.vc.csv) ], one per request
    bivariate_correlations = MPH_REML.out.correlations.map { meta, correlations -> [stripMphPairRouteState(meta), correlations] } // channel: [ val(meta), path(mq.cor.csv) ]
    bivariate_fixed        = MPH_REML.out.fixed_effects.map { meta, fixed_effects -> [stripMphPairRouteState(meta), fixed_effects] } // channel: [ val(meta), path(mq.blue.csv) ]
    bivariate_iterations   = MPH_REML.out.iterations.map { meta, iterations -> [stripMphPairRouteState(meta), iterations] } // channel: [ val(meta), path(mq.iter.csv) ]
    bivariate_log          = MPH_REML.out.log.map { meta, native_log -> [stripMphPairRouteState(meta), native_log] } // channel: [ val(meta), path(log) ]
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
        mph_trait_names: [pair_meta.left_trait_id, pair_meta.right_trait_id],
        mph_effective: buildMphEffectiveSettings(pair_meta.request_options.mph, pair_meta.native_args),
    ]
}

// The staged bundle's prefix, taken from the member whose name defines it. MPH's `--grm_list` names prefixes
// rather than paths, so this is the value the fit is addressed by.
def grmPrefix(grm_files) {
    def grm_iid_name = grm_files.find { grm_file -> grm_file.name.endsWith('.grm.iid') }.name
    return grm_iid_name.substring(0, grm_iid_name.length() - '.grm.iid'.length())
}

// One line of a completed process output read in an operator, which is legal in the way a `path` input read
// inside a `script:` block is not.
def readCsvHeader(csv) {
    return csv.readLines()[0].split(',').toList()
}

// Render the curated options and validated native additions for the configured MPH argument closure.
def buildMphEffectiveSettings(options, native_args) {
    return [
        native_arguments: [
            options.iterations != null ? "--num_iterations ${options.iterations}" : null,
            options.tolerance != null ? "--tolerance ${options.tolerance}" : null,
            options.random_vectors != null ? "--num_random_vectors ${options.random_vectors}" : null,
            options.seed != null ? "--seed ${options.seed}" : null,
            options.save_memory ? '--save_memory' : null,
        ].findAll { argument -> argument != null } + (native_args ?: []).collect { argument -> argument.toString() }
    ]
}

// The route-local keys exist so the atoms, the adapters and the configured publication closures can see the
// selected token and its resolved settings. They are stripped before emission, so a consumer receives the
// request identity it supplied and the route cannot leak its own bookkeeping.
def stripMphPairRouteState(meta) {
    return meta.findAll { name, _value -> !(name in ['matrix_key', 'mph_estimator', 'mph_trait_names', 'mph_effective']) }
}
