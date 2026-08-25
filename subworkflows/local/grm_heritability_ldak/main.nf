// Estimate SNP heritability from a pre-built LDAK kinship matrix.
// Constituent local modules report versions directly to the run-wide topic, so this subworkflow emits no versions.
//

// MODULE: Local to the pipeline
include { LDAK_ADJUSTGRM } from '../../../modules/local/ldak/adjustgrm/main'
include { LDAK_HE        } from '../../../modules/local/ldak/he/main'
include { LDAK_PCGC      } from '../../../modules/local/ldak/pcgc/main'
include { LDAK_REML      } from '../../../modules/local/ldak/reml/main'

workflow GRM_HERITABILITY_LDAK {
    take:
    ch_grm // channel: [ val(meta), path(grm_files) ], pre-built LDAK GRM bundle, once per analysis
    ch_pheno // channel: [ val(meta2), path(phenotype_file), val(prevalence) ], prevalence is [] when unused
    ch_qcovar // channel: [ val(meta3), path(quant_covariates_file) ], use [] for the optional file
    ch_covar // channel: [ val(meta4), path(cat_covariates_file) ], use [] for the optional file
    ch_keep // channel: [ val(meta5), path(keep_file) ], use [] for the optional file
    ch_estimator // channel: [ val(meta6), val(estimator) ], route selector once per analysis
    ch_adjust_covar // channel: [ val(meta7), path(adjustment_covariates_file) ], use [] when absent

    main:
    def ch_estimators = ch_estimator.map { meta, estimator ->
        if (!['reml', 'he', 'pcgc'].contains(estimator)) {
            error("[nf-core/gwas] ERROR: GRM_HERITABILITY_LDAK estimator must be 'reml', 'he' or 'pcgc', got '${estimator}'")
        }
        tuple(meta.id, estimator)
    }

    def ch_analyses = ch_grm
        .map { meta, grm_files -> tuple(meta.id, tuple(meta, grm_files)) }
        .join(
            ch_pheno.map { meta2, phenotype_file, prevalence -> tuple(meta2.id, tuple(meta2, phenotype_file, prevalence)) },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_qcovar.map { meta3, quant_covariates_file -> tuple(meta3.id, tuple(meta3, quant_covariates_file)) },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_covar.map { meta4, cat_covariates_file -> tuple(meta4.id, tuple(meta4, cat_covariates_file)) },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_keep.map { meta5, keep_file -> tuple(meta5.id, tuple(meta5, keep_file)) },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_estimators,
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_adjust_covar.map { meta7, adjustment_covariates_file -> tuple(meta7.id, tuple(meta7, adjustment_covariates_file)) },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )

    // Haseman-Elston and PCGC residualise the phenotype on the covariates but read the kinship matrix
    // as-is, so that matrix must first be regressed on the same covariates. Passing the covariates to
    // both --adjust-grm and --he/--pcgc is LDAK's documented usage, not double-counting: one call
    // projects the kinship, the other residualises the phenotype. REML fits covariates and kinship
    // jointly and therefore needs no adjusted matrix. The keep list travels with the adjustment as
    // well, because LDAK requires the kinship to be regressed on the covariates over exactly the
    // samples the estimator will use rather than over the full cohort.
    def ch_adjust_routes = ch_analyses.branch { _analysis_id, _grm, _pheno, _qcovar, _covar, _keep, estimator, adjustment_covar ->
        adjusted: estimator != 'reml' && adjustment_covar[1]
        direct: true
    }

    def ch_adjustgrm_state = ch_adjust_routes.adjusted.multiMap { analysis_id, grm, pheno, _qcovar, _covar, keep, _estimator, adjustment_covar ->
        grm: tuple([id: "${analysis_id}.adjusted"], grm[1])
        pheno: tuple([id: "${analysis_id}.adjusted"], pheno[1])
        keep: keep
        adjustment_covar: adjustment_covar
        focal_meta: tuple("${analysis_id}.adjusted", grm[0])
    }
    LDAK_ADJUSTGRM(
        ch_adjustgrm_state.grm,
        ch_adjustgrm_state.pheno,
        ch_adjustgrm_state.keep,
        ch_adjustgrm_state.adjustment_covar,
    )

    def ch_adjusted_grm = LDAK_ADJUSTGRM.out.adjusted_grm
        .map { execution_meta, grm_bin, grm_id, grm_details, grm_adjust, grm_root ->
            tuple(execution_meta.id, [grm_bin, grm_id, grm_details, grm_adjust, grm_root])
        }
        .join(
            ch_adjustgrm_state.focal_meta,
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { _execution_id, grm_files, focal_meta -> tuple(focal_meta, grm_files) }

    def ch_adjustgrm_logs = LDAK_ADJUSTGRM.out.log
        .map { execution_meta, log_file -> tuple(execution_meta.id, log_file) }
        .join(
            ch_adjustgrm_state.focal_meta,
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { _execution_id, log_file, focal_meta -> tuple(focal_meta, log_file) }

    def ch_estimator_grm = ch_adjust_routes.direct
        .map { analysis_id, grm, _pheno, _qcovar, _covar, _keep, _estimator, _adjustment_covar ->
            tuple(analysis_id, grm[1], false)
        }
        .mix(ch_adjusted_grm.map { focal_meta, grm_files -> tuple(focal_meta.id, grm_files, true) })

    def ch_invocations = ch_analyses
        .join(
            ch_estimator_grm,
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .branch { _analysis_id, grm, pheno, qcovar, covar, keep, estimator, adjustment_covar, estimator_grm_files, matrix_adjusted ->
            reml: estimator == 'reml'
            return tuple(tuple(grm[0], pheno[1], pheno[2]), tuple(grm[0], estimator_grm_files), keep, qcovar, covar)
            he: estimator == 'he'
            return tuple(
                tuple(grm[0], pheno[1]),
                tuple(grm[0], estimator_grm_files),
                keep,
                matrix_adjusted ? adjustment_covar : qcovar,
                matrix_adjusted ? tuple(covar[0], []) : covar,
            )
            pcgc: estimator == 'pcgc'
            return tuple(
                tuple(grm[0], pheno[1], pheno[2]),
                tuple(grm[0], estimator_grm_files),
                keep,
                matrix_adjusted ? adjustment_covar : qcovar,
                matrix_adjusted ? tuple(covar[0], []) : covar,
            )
        }

    def ch_reml_invocations = ch_invocations.reml.multiMap { pheno, grm, keep, qcovar, covar ->
        pheno: pheno
        grm: grm
        keep: keep
        qcovar: qcovar
        covar: covar
    }
    LDAK_REML(
        ch_reml_invocations.pheno,
        ch_reml_invocations.grm,
        ch_reml_invocations.keep,
        ch_reml_invocations.qcovar,
        ch_reml_invocations.covar,
    )

    def ch_he_invocations = ch_invocations.he.multiMap { pheno, grm, keep, qcovar, covar ->
        pheno: pheno
        grm: grm
        keep: keep
        qcovar: qcovar
        covar: covar
    }
    LDAK_HE(
        ch_he_invocations.pheno,
        ch_he_invocations.grm,
        ch_he_invocations.keep,
        ch_he_invocations.qcovar,
        ch_he_invocations.covar,
    )

    def ch_pcgc_invocations = ch_invocations.pcgc.multiMap { pheno, grm, keep, qcovar, covar ->
        pheno: pheno
        grm: grm
        keep: keep
        qcovar: qcovar
        covar: covar
    }
    LDAK_PCGC(
        ch_pcgc_invocations.pheno,
        ch_pcgc_invocations.grm,
        ch_pcgc_invocations.keep,
        ch_pcgc_invocations.qcovar,
        ch_pcgc_invocations.covar,
    )

    def ch_coeff = LDAK_REML.out.coeff
    ch_coeff = ch_coeff.mix(LDAK_HE.out.coeff)
    ch_coeff = ch_coeff.mix(LDAK_PCGC.out.coeff)

    def ch_cross = LDAK_REML.out.cross
    ch_cross = ch_cross.mix(LDAK_HE.out.cross)
    ch_cross = ch_cross.mix(LDAK_PCGC.out.cross)

    def ch_share = LDAK_REML.out.share
    ch_share = ch_share.mix(LDAK_HE.out.share)
    ch_share = ch_share.mix(LDAK_PCGC.out.share)

    def ch_progress = LDAK_REML.out.progress
    ch_progress = ch_progress.mix(LDAK_HE.out.progress)
    ch_progress = ch_progress.mix(LDAK_PCGC.out.progress)

    def ch_logs = LDAK_REML.out.log
    ch_logs = ch_logs.mix(LDAK_HE.out.log)
    ch_logs = ch_logs.mix(LDAK_PCGC.out.log)
    ch_logs = ch_logs.mix(ch_adjustgrm_logs)

    emit:
    reml_results   = LDAK_REML.out.reml_results // channel: [ val(meta), path(reml) ], REML-route records only
    reml_liability = LDAK_REML.out.reml_liability // channel: [ val(meta), path(reml_liab) ], REML records given a prevalence
    he_results     = LDAK_HE.out.he_results // channel: [ val(meta), path(he) ], Haseman-Elston-route records only
    pcgc_results   = LDAK_PCGC.out.pcgc_results // channel: [ val(meta), path(pcgc) ], PCGC-route records only
    pcgc_marginal  = LDAK_PCGC.out.pcgc_marginal // channel: [ val(meta), path(pcgc_marginal) ], optional PCGC-route records
    adjusted_grm   = ch_adjusted_grm // channel: [ val(meta), path(grm_files) ], adjusted-route records only
    coeff          = ch_coeff // channel: [ val(meta), path(coeff) ]
    cross          = ch_cross // channel: [ val(meta), path(cross) ]
    share          = ch_share // channel: [ val(meta), path(share) ]
    progress       = ch_progress // channel: [ val(meta), path(progress) ]
    logs           = ch_logs // channel: [ val(meta), path(log) ]
}
