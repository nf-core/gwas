// Estimate heritability with GCTA GREML or GREML-LDMS from prepared relatedness matrices.
// Both processes report on the run-wide versions topic, so this subworkflow emits no versions.

// MODULE: Local to the pipeline
include { GCTA_REML     } from '../../../modules/local/gcta/reml/main'
include { GCTA_REMLLDMS } from '../../../modules/local/gcta/remlldms/main'

workflow GRM_HERITABILITY_GCTA {
    take:
    ch_grm // channel: [ val(meta), path(grm_files), val(grm_prefixes) ], prefixes are [] for greml
    ch_pheno // channel: [ val(meta2), path(phenotypes_file) ], once per analysis
    ch_qcovar // channel: [ val(meta3), path(quant_covariates_file) ], use [] for the optional file
    ch_covar // channel: [ val(meta4), path(cat_covariates_file) ], use [] for the optional file
    ch_estimator // channel: [ val(meta5), val(estimator) ], 'greml' or 'greml_ldms'

    main:
    ch_estimators = ch_estimator.map { meta5, estimator ->
        if (!['greml', 'greml_ldms'].contains(estimator)) {
            error("[nf-core/gwas] ERROR: GRM_HERITABILITY_GCTA estimator must be 'greml' or 'greml_ldms', got '${estimator}'")
        }
        tuple([meta5.id, meta5.gcta_estimator], estimator)
    }

    ch_routes = ch_grm
        .map { meta, grm_files, grm_prefixes -> tuple([meta.id, meta.gcta_estimator], tuple(meta, grm_files, grm_prefixes)) }
        .join(
            ch_pheno.map { meta2, phenotypes_file -> tuple([meta2.id, meta2.gcta_estimator], tuple(meta2, phenotypes_file)) },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_qcovar.map { meta3, quant_covariates_file -> tuple([meta3.id, meta3.gcta_estimator], tuple(meta3, quant_covariates_file)) },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_covar.map { meta4, cat_covariates_file -> tuple([meta4.id, meta4.gcta_estimator], tuple(meta4, cat_covariates_file)) },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(ch_estimators, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { route_key, grm, pheno, qcovar, covar, estimator -> tuple(route_key, grm, pheno, qcovar, covar, estimator) }
        .branch { _route_key, grm, pheno, qcovar, covar, estimator ->
            greml: estimator == 'greml'
            return tuple(tuple(grm[0], grm[1]), pheno, qcovar, covar)
            greml_ldms: estimator == 'greml_ldms'
            return tuple(grm, pheno, qcovar, covar)
        }

    ch_greml_inputs = ch_routes.greml.multiMap { grm, pheno, qcovar, covar ->
        grm: grm
        pheno: pheno
        qcovar: qcovar
        covar: covar
    }
    GCTA_REML(ch_greml_inputs.grm, ch_greml_inputs.pheno, ch_greml_inputs.qcovar, ch_greml_inputs.covar)

    ch_greml_ldms_inputs = ch_routes.greml_ldms.multiMap { grm, pheno, qcovar, covar ->
        grm: grm
        pheno: pheno
        qcovar: qcovar
        covar: covar
    }
    GCTA_REMLLDMS(ch_greml_ldms_inputs.grm, ch_greml_ldms_inputs.pheno, ch_greml_ldms_inputs.qcovar, ch_greml_ldms_inputs.covar)

    ch_heritability = GCTA_REML.out.reml_results.mix(GCTA_REMLLDMS.out.reml_results)
    ch_logs = GCTA_REML.out.log.mix(GCTA_REMLLDMS.out.log)

    emit:
    heritability = ch_heritability // channel: [ val(meta), path(hsq) ], one record per analysis
    logs         = ch_logs // channel: [ val(meta), path(log) ]
}
