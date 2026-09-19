// Route unary analysis units through GCTA GREML/GREML-LDMS and LDAK REML/HE/PCGC.
// Relatedness construction remains on the pipeline spine; this controller owns estimator adaptation and
// covariate-adjusted LDAK derivative reuse above the estimator fan-out.

// MODULE: Local to the pipeline
include { GCTA_REML                  } from '../../../modules/local/gcta/reml/main'
include { GCTA_REMLLDMS              } from '../../../modules/local/gcta/remlldms/main'
include { LDAK_ADJUSTGRM             } from '../../../modules/local/ldak/adjustgrm/main'
include { LDAK_REML                  } from '../../../modules/local/ldak/reml/main'
include { LDAK_HE                    } from '../../../modules/local/ldak/he/main'
include { LDAK_PCGC                  } from '../../../modules/local/ldak/pcgc/main'

// FUNCTION: Local to the pipeline
include { buildScientificArtifactKey } from '../utils_nfcore_gwas_pipeline'
include { digestFileBytes            } from '../utils_nfcore_gwas_pipeline'
include { digestIdentityText         } from '../utils_nfcore_gwas_pipeline'

workflow ROUTE_GRM_HERITABILITY {
    take:
    ch_dense_matrices // channel: [ val(meta), [ path(grm_file), ... ] ], dense GCTA matrices fanned to unary consumers
    ch_ldms_matrices // channel: [ val(meta), [ path(grm_file), ... ], val(grm_prefixes) ], LDMS matrix families fanned to unary consumers
    ch_ldak_kinship_matrices // channel: [ val(meta), val(artifact_key), [ path(grm_file), ... ], path(keep) ], selected LDAK base/subset artifact
    ch_headerless_phenotypes // channel: [ val(meta), path(phenotype), path(quant_covariates), path(cat_covariates) ], optional covariates are []
    ch_adjustment_covariates // channel: [ val(meta), path(adjustment_covariates) ], only analyses with covariates

    main:

    // GCTA GREML and GREML-LDMS preserve their native matrix families and result contracts.
    def ch_greml_matrices = ch_dense_matrices
        .map { meta, grm_files -> [meta, grm_files, [], 'greml'] }
        .mix(
            ch_ldms_matrices.map { meta, grm_files, grm_prefixes -> [meta, grm_files, grm_prefixes, 'greml_ldms'] }
        )

    def ch_greml_routes = ch_greml_matrices
        .combine(ch_headerless_phenotypes, by: 0)
        .branch { _meta, _grm_files, _grm_prefixes, estimator, _phenotype, _quant_covariates, _cat_covariates ->
            greml: estimator == 'greml'
            greml_ldms: estimator == 'greml_ldms'
        }

    def ch_greml_inputs = ch_greml_routes.greml.multiMap { meta, grm_files, _grm_prefixes, estimator, phenotype, quant_covariates, cat_covariates ->
        def route_meta = meta + [gcta_estimator: estimator]
        grm: [route_meta, grm_files]
        pheno: [route_meta, phenotype]
        qcovar: [route_meta, quant_covariates]
        covar: [route_meta, cat_covariates]
    }
    GCTA_REML(ch_greml_inputs.grm, ch_greml_inputs.pheno, ch_greml_inputs.qcovar, ch_greml_inputs.covar)

    // GCTA's component order is scientific result identity; order complete staged bundles by the planner's
    // declared family before the native module derives its MGRM entries from their filenames. The declared
    // family and the staged bundle are two members of one `PLINK_PREPARE_GRM_LDMS_GCTA.out.grm_family`
    // emission, so every declared prefix has its files here and the ordering never drops a component.
    def ch_ldms_inputs = ch_greml_routes.greml_ldms.multiMap { meta, grm_files, grm_prefixes, estimator, phenotype, quant_covariates, cat_covariates ->
        def route_meta = meta + [gcta_estimator: estimator]
        grm: [route_meta, grm_prefixes.collectMany { grm_prefix -> grm_files.findAll { grm_file -> grm_file.name.startsWith("${grm_prefix}.grm.") } }]
        pheno: [route_meta, phenotype]
        qcovar: [route_meta, quant_covariates]
        covar: [route_meta, cat_covariates]
    }
    GCTA_REMLLDMS(ch_ldms_inputs.grm, ch_ldms_inputs.pheno, ch_ldms_inputs.qcovar, ch_ldms_inputs.covar)

    // Expand one selected LDAK parent per analysis into its requested estimators while retaining the parent
    // artifact key beside, rather than inside, focal metadata.
    // The remainder join also emits covariate-only records for other methods; their absent parent is filtered.
    def ch_ldak_inputs = ch_ldak_kinship_matrices
        .join(ch_headerless_phenotypes, failOnDuplicate: true)
        .join(ch_adjustment_covariates, remainder: true)
        .filter { record -> record[1] != null }
        .flatMap { meta, parent_key, grm_files, keep, phenotype, quant_covariates, cat_covariates, adjustment_covariates ->
            [
                [method: 'ldak_reml', estimator: 'reml'],
                [method: 'ldak_he', estimator: 'he'],
                [method: 'ldak_pcgc', estimator: 'pcgc'],
            ].findAll { route -> route.method in meta.heritability_methods }.collect { route ->
                [meta, parent_key, grm_files, keep, phenotype, quant_covariates, cat_covariates, adjustment_covariates ?: [], route.estimator]
            }
        }

    // Adjustment identity is independent of the estimator. HE and PCGC consumers with the same parent,
    // effective sample subset, numerical design, and native options therefore reduce to one native task.
    def ch_adjustment_consumers = ch_ldak_inputs
        .filter { _meta, _parent_key, _grm_files, _keep, _phenotype, _quant_covariates, _cat_covariates, adjustment_covariates, estimator ->
            estimator in ['he', 'pcgc'] && adjustment_covariates
        }
        .map { meta, parent_key, grm_files, keep, phenotype, quant_covariates, cat_covariates, adjustment_covariates, estimator ->
            def request = buildLdakAdjustmentRequest(parent_key, phenotype, keep, adjustment_covariates)
            [request.key, meta, request, grm_files, keep, phenotype, quant_covariates, cat_covariates, adjustment_covariates, estimator]
        }

    def ch_adjustment_builds = ch_adjustment_consumers
        .map { key, _meta, request, grm_files, keep, phenotype, _quant_covariates, _cat_covariates, adjustment_covariates, _estimator ->
            [key, request, grm_files, keep, phenotype, adjustment_covariates]
        }
        .unique { key, _request, _grm_files, _keep, _phenotype, _adjustment_covariates -> key }
        .multiMap { _key, request, grm_files, keep, phenotype, adjustment_covariates ->
            def artifact_meta = buildLdakAdjustmentArtifactMeta(request)
            grm: [artifact_meta, grm_files]
            pheno: [artifact_meta, phenotype]
            keep: [artifact_meta, keep]
            adjustment_covariates: [artifact_meta, adjustment_covariates]
        }

    // LDAK records the covariate filename in the adjusted artifact's native `.root` contract. Every consumer
    // of one byte-identical adjustment artifact must therefore stage the same representative input file that
    // constructed it, rather than a focal-analysis copy with another basename.
    // .grm.root records the --covar path as typed and
    // --he/--pcgc exit 1 on any other path (LDAK 6 and 6.3). Retire if LDAK stops checking the root
    // path or --check-root NO is adopted as policy.
    def ch_adjustment_reference_covariates = ch_adjustment_consumers
        .map { key, _meta, _request, _grm_files, _keep, _phenotype, _quant_covariates, _cat_covariates, adjustment_covariates, _estimator ->
            [key, adjustment_covariates]
        }
        .unique { key, _adjustment_covariates -> key }

    LDAK_ADJUSTGRM(
        ch_adjustment_builds.grm,
        ch_adjustment_builds.pheno,
        ch_adjustment_builds.keep,
        ch_adjustment_builds.adjustment_covariates,
    )

    def ch_adjusted_grm = LDAK_ADJUSTGRM.out.adjusted_grm
        .map { artifact_meta, grm_files ->
            [artifact_meta.key, grm_files]
        }
        .combine(ch_adjustment_reference_covariates, by: 0)

    def ch_adjusted_invocations = ch_adjustment_consumers
        .map { key, meta, _request, _grm_files, keep, phenotype, _quant_covariates, _cat_covariates, adjustment_covariates, estimator ->
            [key, meta, keep, phenotype, adjustment_covariates, estimator]
        }
        .combine(ch_adjusted_grm, by: 0)
        .map { _key, meta, keep, phenotype, _adjustment_covariates, estimator, adjusted_grm_files, artifact_covariates ->
            [meta, estimator, adjusted_grm_files, keep, phenotype, artifact_covariates, []]
        }

    // REML always uses the selected unadjusted parent. HE/PCGC do likewise only when there is no covariate
    // design requiring an adjusted child.
    def ch_direct_invocations = ch_ldak_inputs
        .filter { _meta, _parent_key, _grm_files, _keep, _phenotype, _quant_covariates, _cat_covariates, adjustment_covariates, estimator ->
            estimator == 'reml' || !adjustment_covariates
        }
        .map { meta, _parent_key, grm_files, keep, phenotype, quant_covariates, cat_covariates, _adjustment_covariates, estimator ->
            [meta, estimator, grm_files, keep, phenotype, quant_covariates, cat_covariates]
        }

    def ch_ldak_invocations = ch_direct_invocations
        .mix(ch_adjusted_invocations)
        .branch { _meta, estimator, _grm_files, _keep, _phenotype, _quant_covariates, _cat_covariates ->
            reml: estimator == 'reml'
            he: estimator == 'he'
            pcgc: estimator == 'pcgc'
        }

    def ch_ldak_reml_inputs = ch_ldak_invocations.reml.multiMap { meta, _estimator, grm_files, keep, phenotype, quant_covariates, cat_covariates ->
        pheno: [meta, phenotype]
        grm: [meta, grm_files]
        keep: [meta, keep ?: []]
        qcovar: [meta, quant_covariates]
        covar: [meta, cat_covariates]
    }
    LDAK_REML(ch_ldak_reml_inputs.pheno, ch_ldak_reml_inputs.grm, ch_ldak_reml_inputs.keep, ch_ldak_reml_inputs.qcovar, ch_ldak_reml_inputs.covar)

    def ch_ldak_he_inputs = ch_ldak_invocations.he.multiMap { meta, _estimator, grm_files, keep, phenotype, quant_covariates, cat_covariates ->
        pheno: [meta, phenotype]
        grm: [meta, grm_files]
        keep: [meta, keep ?: []]
        qcovar: [meta, quant_covariates]
        covar: [meta, cat_covariates]
    }
    LDAK_HE(ch_ldak_he_inputs.pheno, ch_ldak_he_inputs.grm, ch_ldak_he_inputs.keep, ch_ldak_he_inputs.qcovar, ch_ldak_he_inputs.covar)

    def ch_ldak_pcgc_inputs = ch_ldak_invocations.pcgc.multiMap { meta, _estimator, grm_files, keep, phenotype, quant_covariates, cat_covariates ->
        pheno: [meta, phenotype, meta.population_prevalence]
        grm: [meta, grm_files]
        keep: [meta, keep ?: []]
        qcovar: [meta, quant_covariates]
        covar: [meta, cat_covariates]
    }
    LDAK_PCGC(ch_ldak_pcgc_inputs.pheno, ch_ldak_pcgc_inputs.grm, ch_ldak_pcgc_inputs.keep, ch_ldak_pcgc_inputs.qcovar, ch_ldak_pcgc_inputs.covar)

    emit:
    gcta_heritability   = GCTA_REML.out.reml_results.mix(GCTA_REMLLDMS.out.reml_results) // channel: [ val(meta), path(hsq) ], one per analysis per selected GCTA estimator
    ldak_reml_results   = LDAK_REML.out.reml_results // channel: [ val(meta), path(reml) ], one per analysis selecting ldak_reml
    ldak_reml_liability = LDAK_REML.out.reml_liability // channel: [ val(meta), path(reml_liab) ], only for a row declaring population prevalence
    ldak_he_results     = LDAK_HE.out.he_results // channel: [ val(meta), path(he) ], one per analysis selecting ldak_he
    ldak_pcgc_results   = LDAK_PCGC.out.pcgc_results // channel: [ val(meta), path(pcgc) ], one per analysis selecting ldak_pcgc
    ldak_pcgc_marginal  = LDAK_PCGC.out.pcgc_marginal // channel: [ val(meta), path(pcgc_marginal) ], optional PCGC-route records
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def getLdakAdjustmentSampleIdentity(phenotype_file, keep_file) {
    def phenotype_path = phenotype_file instanceof java.nio.file.Path ? phenotype_file : phenotype_file.toPath()
    def phenotype_samples = java.nio.file.Files.readAllLines(phenotype_path).collect { line -> line.split('\t', -1) }.findAll { fields -> fields[2] != 'NA' }.collect { fields -> "${fields[0]}\t${fields[1]}" } as Set
    def effective_samples = keep_file
        ? java.nio.file.Files.readAllLines(keep_file instanceof java.nio.file.Path ? keep_file : keep_file.toPath()).collect { line -> line.split(/\s+/, -1).take(2).join('\t') }.findAll { sample -> sample in phenotype_samples }
        : phenotype_samples.toList()
    return digestIdentityText(effective_samples.sort().join('\n'))
}

def buildLdakAdjustmentRequest(parent_key, phenotype_file, keep_file, adjustment_covariates_file) {
    def settings = [
        sample_subset: getLdakAdjustmentSampleIdentity(phenotype_file, keep_file),
        adjustment_covariates: [sha256: digestFileBytes(adjustment_covariates_file)],
        native_options: [],
    ]
    return [
        layer: 'derived',
        type: 'ldak_adjusted',
        parent_key: parent_key,
        settings: settings,
        key: buildScientificArtifactKey([layer: 'derived', type: 'ldak_adjusted', parent_key: parent_key], settings),
    ]
}

def buildLdakAdjustmentArtifactMeta(request) {
    return [
        id: "${request.type}.${request.key}",
        key: request.key,
        layer: request.layer,
        artifact_type: request.type,
        parent_key: request.parent_key,
        settings: request.settings,
    ]
}
