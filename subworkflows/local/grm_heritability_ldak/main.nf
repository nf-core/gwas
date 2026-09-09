// Estimate SNP heritability from a selected pre-built LDAK kinship artifact.
// Constituent local modules report versions directly to the run-wide topic, so this subworkflow emits no versions.

// MODULE: Local to the pipeline
include { LDAK_HE   } from '../../../modules/local/ldak/he/main'
include { LDAK_PCGC } from '../../../modules/local/ldak/pcgc/main'
include { LDAK_REML } from '../../../modules/local/ldak/reml/main'

workflow GRM_HERITABILITY_LDAK {
    take:
    ch_grm // channel: [ val(meta), path(grm_files) ], selected base or adjusted LDAK GRM
    ch_pheno // channel: [ val(meta2), path(phenotype_file), val(prevalence) ], prevalence is [] when unused
    ch_qcovar // channel: [ val(meta3), path(quant_covariates_file) ], use [] for the optional file
    ch_covar // channel: [ val(meta4), path(cat_covariates_file) ], use [] for the optional file
    ch_keep // channel: [ val(meta5), path(keep_file) ], use [] for the optional file
    ch_estimator // channel: [ val(meta6), val(estimator) ], route selector once per analysis

    main:
    def ch_estimators = ch_estimator.map { meta, estimator ->
        if (!['reml', 'he', 'pcgc'].contains(estimator)) {
            error("[nf-core/gwas] ERROR: GRM_HERITABILITY_LDAK estimator must be 'reml', 'he' or 'pcgc', got '${estimator}'")
        }
        [meta.id, estimator]
    }

    def ch_analyses = ch_grm
        .map { meta, grm_files -> [meta.id, [meta, grm_files]] }
        .join(
            ch_pheno.map { meta2, phenotype_file, prevalence -> [meta2.id, [meta2, phenotype_file, prevalence]] },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_qcovar.map { meta3, quant_covariates_file -> [meta3.id, [meta3, quant_covariates_file]] },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_covar.map { meta4, cat_covariates_file -> [meta4.id, [meta4, cat_covariates_file]] },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_keep.map { meta5, keep_file -> [meta5.id, [meta5, keep_file]] },
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

    def ch_invocations = ch_analyses.branch { _analysis_id, grm, pheno, qcovar, covar, keep, estimator ->
        reml: estimator == 'reml'
        return [[grm[0], pheno[1], pheno[2]], [grm[0], grm[1]], keep, qcovar, covar]
        he: estimator == 'he'
        return [[grm[0], pheno[1], pheno[2]], [grm[0], grm[1]], keep, qcovar, covar]
        pcgc: estimator == 'pcgc'
        return [[grm[0], pheno[1], pheno[2]], [grm[0], grm[1]], keep, qcovar, covar]
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

    emit:
    reml_results   = LDAK_REML.out.reml_results // channel: [ val(meta), path(reml) ], REML-route records only
    reml_liability = LDAK_REML.out.reml_liability // channel: [ val(meta), path(reml_liab) ], REML records given a prevalence
    he_results     = LDAK_HE.out.he_results // channel: [ val(meta), path(he) ], Haseman-Elston-route records only
    pcgc_results   = LDAK_PCGC.out.pcgc_results // channel: [ val(meta), path(pcgc) ], PCGC-route records only
    pcgc_marginal  = LDAK_PCGC.out.pcgc_marginal // channel: [ val(meta), path(pcgc_marginal) ], optional PCGC-route records
    coeff          = ch_coeff // channel: [ val(meta), path(coeff) ]
    cross          = ch_cross // channel: [ val(meta), path(cross) ]
    share          = ch_share // channel: [ val(meta), path(share) ]
    progress       = ch_progress // channel: [ val(meta), path(progress) ]
    logs           = ch_logs // channel: [ val(meta), path(log) ]
}
