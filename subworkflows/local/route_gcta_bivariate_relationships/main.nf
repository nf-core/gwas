// Route nf-core/gwas individual-level relationship requests through the GCTA bivariate estimators. A
// relationship is a declared, oriented pair of analysis units over one cohort; it is a different domain from a
// summary-statistics pair request and the two are never merged here. Dense REML, REML-LDMS, dense HEreg and
// HEreg-LDMS share one controller because they share the relationship definition and the one ordered two-trait
// table prepared from it on the spine: a relationship
// selecting several of them is prepared exactly once, and the dense or LDMS matrix a REML request already
// caused to be built is the same matrix its HEreg sibling consumes, because the matrix reuse key is derived
// only from the cohort, the genotype bundle and the declared matrix settings.
//
// This is pipeline routing, reuse identity and scientific/publication policy, not an nf-core/modules submission
// candidate. Declared left/right orientation is the execution and presentation order and is taken as given: a
// reversed duplicate is refused upstream by validation, and nothing here re-derives, reorders or de-duplicates
// an orientation. That contract is load-bearing for HEreg in particular, whose cross-product coefficient is
// fitted on the lower triangle only and therefore changes when the two traits swap sides. Every constituent
// process reports directly to the run-wide versions topic, so this subworkflow emits no versions, and it reads
// no params, no workflow and no projectDir.

// MODULES: Installed directly from nf-core/modules
include { GCTA_BIVARIATEREML      } from '../../../modules/nf-core/gcta/bivariatereml/main'
include { GCTA_BIVARIATEREMLLDMS  } from '../../../modules/nf-core/gcta/bivariateremlldms/main'

// MODULE: Local to the pipeline
include { GCTA_BIVARIATEHEREG     } from '../../../modules/local/gcta/bivariatehereg/main'
include { GCTA_BIVARIATEHEREGLDMS } from '../../../modules/local/gcta/bivariateheregldms/main'

workflow ROUTE_GCTA_BIVARIATE_RELATIONSHIPS {
    take:
    ch_relationships // channel: [ val(meta), [ path(genotype_file), ... ], path(pair_quant_covariates), path(pair_cat_covariates) ], one validated row per selected method per relationship
    ch_prepared_relationship_traits // channel: [ val(relationship_meta), path(phenotype), path(quant_covariates), path(cat_covariates) ], headerless, one per relationship; [] for an absent covariate table
    ch_relationship_dense_matrices // channel: [ val(meta), [ path(grm_file), ... ] ], the dense GCTA matrices built for relationship requests, already fanned out per request
    ch_relationship_ldms_matrices // channel: [ val(meta), [ path(grm_file), ... ], val(grm_prefixes) ], ordered LDMS matrix families built for relationship requests

    main:

    // The ordered two-trait table and the normalised pair covariates are the scientific pair, prepared once
    // per relationship on the spine because the MPH pair controller consumes the same artifact. This
    // controller fans that one record back out by relationship ID, so selecting several GCTA estimators still
    // costs one preparation and defines no second endpoint sample set.
    def ch_prepared_pairs = ch_relationships
        .map { meta, _genotype_files, _pair_quant_covariates, _pair_cat_covariates -> [meta.relationship_id, meta] }
        .combine(
            ch_prepared_relationship_traits.map { relationship_meta, phenotype, quant_covariates, cat_covariates ->
                [relationship_meta.relationship_id, phenotype, quant_covariates ?: [], cat_covariates ?: []]
            },
            by: 0,
        )
        .map { _relationship_id, meta, phenotype, quant_covariates, cat_covariates ->
            [meta.request_id, meta, phenotype, quant_covariates, cat_covariates]
        }

    //
    // Dense GCTA relationship requests: bivariate REML and bivariate HEreg
    //
    // The request identity remains the primary metadata ID while the matrix content key records reuse
    // attribution. Both dense estimators are addressed from the same prepared pair and the same matrix stream,
    // so the branch below is the only place they diverge.
    def ch_bivariate_matrices = ch_relationship_dense_matrices
        .filter { meta, _grm_files -> meta.matrix_kind == 'gcta_dense' }
        .map { meta, grm_files ->
            [meta.request_id, meta, grm_files]
        }

    def ch_dense_prepared_pairs = ch_prepared_pairs.filter { _request_id, pair_meta, _phenotype, _quant_covariates, _cat_covariates -> pair_meta.matrix_kind == 'gcta_dense' }

    def ch_dense_requests = ch_bivariate_matrices
        .join(ch_dense_prepared_pairs, failOnDuplicate: true, failOnMismatch: true)
        .map { _request_id, matrix_meta, grm_files, pair_meta, phenotype, quant_covariates, cat_covariates ->
            [resolveRouteMeta(matrix_meta, pair_meta), grm_files, phenotype, quant_covariates, cat_covariates]
        }
        .branch { route_meta, _grm_files, _phenotype, _quant_covariates, _cat_covariates ->
            reml: route_meta.method == 'gcta_bivariate_reml'
            hereg: route_meta.method == 'gcta_bivariate_he'
        }

    //
    // MODULE: primary dense GCTA bivariate REML relationship request
    //
    def ch_bivariate_invocations = ch_dense_requests.reml.multiMap { route_meta, grm_files, phenotype, quant_covariates, cat_covariates ->
        grm: [route_meta, grm_files]
        pheno: [route_meta, phenotype, 1, 2]
        qcovar: [route_meta, quant_covariates]
        covar: [route_meta, cat_covariates]
    }

    GCTA_BIVARIATEREML(
        ch_bivariate_invocations.grm,
        ch_bivariate_invocations.pheno,
        ch_bivariate_invocations.qcovar,
        ch_bivariate_invocations.covar,
    )

    //
    // MODULE: primary dense GCTA bivariate HEreg relationship request
    //
    // `gcta --HEreg-bivar` has no covariate parameter: GCTA 1.94.1 lists `--qcovar` and `--covar` among its
    // accepted options and then never reads them, so passing the prepared covariate tables would silently
    // return an unadjusted estimate. Validation refuses a covariate-bearing HE request before execution and
    // the prepared covariate tables are deliberately dropped here rather than forwarded.
    def ch_hereg_invocations = ch_dense_requests.hereg.multiMap { route_meta, grm_files, phenotype, _quant_covariates, _cat_covariates ->
        grm: [route_meta, grm_files]
        pheno: [route_meta, phenotype, 1, 2]
    }

    GCTA_BIVARIATEHEREG(
        ch_hereg_invocations.grm,
        ch_hereg_invocations.pheno,
    )

    //
    // LDMS GCTA relationship requests: bivariate REML-LDMS and bivariate HEreg-LDMS
    //
    // The request ID and matrix content key remain separate attribution fields so a unary GREML-LDMS request
    // and a pair request can share one scientifically identical ordered matrix family without sharing result
    // identity. The consuming task writes its MGRM control list from the supplied prefix order.
    def ch_bivariate_ldms_matrices = ch_relationship_ldms_matrices
        .filter { meta, _grm_files, _grm_prefixes -> meta.matrix_kind == 'gcta_ldms' }
        .map { meta, grm_files, grm_prefixes ->
            [meta.request_id, meta, grm_files, grm_prefixes]
        }

    def ch_ldms_prepared_pairs = ch_prepared_pairs.filter { _request_id, pair_meta, _phenotype, _quant_covariates, _cat_covariates -> pair_meta.matrix_kind == 'gcta_ldms' }

    def ch_ldms_requests = ch_bivariate_ldms_matrices
        .join(ch_ldms_prepared_pairs, failOnDuplicate: true, failOnMismatch: true)
        .map { _request_id, matrix_meta, grm_files, grm_prefixes, pair_meta, phenotype, quant_covariates, cat_covariates ->
            [resolveLdmsRouteMeta(matrix_meta, pair_meta), grm_files, grm_prefixes, phenotype, quant_covariates, cat_covariates]
        }
        .branch { route_meta, _grm_files, _grm_prefixes, _phenotype, _quant_covariates, _cat_covariates ->
            reml: route_meta.method == 'gcta_bivariate_reml_ldms'
            hereg: route_meta.method == 'gcta_bivariate_he_ldms'
        }

    //
    // MODULE: primary GCTA bivariate REML-LDMS relationship request
    //
    def ch_bivariate_ldms_invocations = ch_ldms_requests.reml.multiMap { route_meta, grm_files, grm_prefixes, phenotype, quant_covariates, cat_covariates ->
        mgrm: [route_meta, grm_files, grm_prefixes]
        pheno: [route_meta, phenotype, 1, 2]
        qcovar: [route_meta, quant_covariates]
        covar: [route_meta, cat_covariates]
    }

    GCTA_BIVARIATEREMLLDMS(
        ch_bivariate_ldms_invocations.mgrm,
        ch_bivariate_ldms_invocations.pheno,
        ch_bivariate_ldms_invocations.qcovar,
        ch_bivariate_ldms_invocations.covar,
    )

    //
    // MODULE: primary GCTA bivariate HEreg-LDMS relationship request
    //
    // The multi-component moment fit emits its own native component and total results.
    def ch_hereg_ldms_invocations = ch_ldms_requests.hereg.multiMap { route_meta, grm_files, grm_prefixes, phenotype, _quant_covariates, _cat_covariates ->
        mgrm: [route_meta, grm_files, grm_prefixes]
        pheno: [route_meta, phenotype, 1, 2]
    }

    GCTA_BIVARIATEHEREGLDMS(
        ch_hereg_ldms_invocations.mgrm,
        ch_hereg_ldms_invocations.pheno,
    )

    emit:
    reml_results       = GCTA_BIVARIATEREML.out.bivariate_results.map { meta, result -> [stripMatrixReuseIdentity(meta), result] } // channel: [ val(meta), path(native.hsq) ]
    reml_log           = GCTA_BIVARIATEREML.out.log_file.map { meta, log -> [stripMatrixReuseIdentity(meta), log] } // channel: [ val(meta), path(native.log) ]
    hereg_results      = GCTA_BIVARIATEHEREG.out.hereg_results.map { meta, result -> [stripMatrixReuseIdentity(meta), result] } // channel: [ val(meta), path(native.HEreg) ]
    hereg_log          = GCTA_BIVARIATEHEREG.out.log.map { meta, log -> [stripMatrixReuseIdentity(meta), log] } // channel: [ val(meta), path(native.log) ]
    reml_ldms_results  = GCTA_BIVARIATEREMLLDMS.out.bivariate_results.map { meta, result -> [stripNativeMatrixIdentity(meta), result] } // channel: [ val(meta), path(native.hsq) ]
    reml_ldms_log      = GCTA_BIVARIATEREMLLDMS.out.log_file.map { meta, log -> [stripNativeMatrixIdentity(meta), log] } // channel: [ val(meta), path(native.log) ]
    hereg_ldms_results = GCTA_BIVARIATEHEREGLDMS.out.hereg_results.map { meta, result -> [stripNativeMatrixIdentity(meta), result] } // channel: [ val(meta), path(native.HEreg) ]
    hereg_ldms_log     = GCTA_BIVARIATEHEREGLDMS.out.log.map { meta, log -> [stripNativeMatrixIdentity(meta), log] } // channel: [ val(meta), path(native.log) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Fold the matrix reuse identity into the request record. Dense atoms resolve the native GRM basename from the
// staged bundle, so the request identity remains primary throughout execution.
def resolveRouteMeta(matrix_meta, pair_meta) {
    return pair_meta + [matrix_key: matrix_meta.matrix_key]
}

def resolveLdmsRouteMeta(matrix_meta, pair_meta) {
    return pair_meta + [id: matrix_meta.matrix_key, matrix_key: matrix_meta.matrix_key]
}

// Dense inputs carry the matrix reuse key for attribution, but it is not part of the emitted request record.
def stripMatrixReuseIdentity(meta) {
    return meta.findAll { key, _value -> key != 'matrix_key' }
}

// LDMS execution uses the matrix reuse key as its task-local output identity. Drop that key and restore the
// request identity before emitting the result.
def stripNativeMatrixIdentity(meta) {
    return stripMatrixReuseIdentity(meta) + [id: meta.request_id]
}
