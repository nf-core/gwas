// Route unary analysis units through LDAK's direct-genotype randomised Haseman-Elston and PCGC estimators.
// No relatedness matrix is requested or built: the estimators read the cohort's PLINK 1 derivative, project
// covariates internally, and take the model, power and weights policy from the same LDAK method options the
// kinship routes use. Emits no versions; reads no params.

// MODULE: Local to the pipeline
include { LDAK_FASTHE                     } from '../../../modules/local/ldak/fasthe/main'
include { LDAK_FASTPCGC                   } from '../../../modules/local/ldak/fastpcgc/main'

// FUNCTION: Local to the pipeline
include { getMethodCapability             } from '../validate_gwas_input/method_registry'
include { getMethodTokensWithCapabilities } from '../validate_gwas_input/method_registry'

workflow ROUTE_LDAK_DIRECT_HERITABILITY {
    take:
    ch_plink1_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam) ], unary analyses whose selectors need PLINK 1
    ch_headerless_phenotypes // channel: [ val(meta), path(phenotype), path(quant_covariates), path(cat_covariates) ], optional covariates are []
    ch_ldak_weights // channel: [ val(meta), path(weights) ], [] unless ldak.weights_policy is 'provided'

    main:

    // A capability query, not a token list: any LDAK heritability estimator that reads the direct PLINK 1
    // stream is routed here, and the native mode it runs in is derived from its declared estimator family.
    def direct_methods = getMethodTokensWithCapabilities([domain: 'heritability', option_family: 'ldak', input_backend: 'direct_plink1_genotypes'])
    def native_modes = [moment_he: 'he', pcgc: 'pcgc']
    def selects_direct = { meta -> meta.heritability_methods.any { method -> method in direct_methods } }

    // Every stream is narrowed by the same predicate before the guarded joins, so a phenotype or weights
    // record belonging to an analysis this route does not serve is not a mismatch. The join key is the
    // analysis identifier, never the meta map, because the three streams carry the same map from three
    // different producers.
    def ch_direct_inputs = ch_plink1_genotypes
        .filter { meta, _bed, _bim, _fam -> selects_direct.call(meta) }
        .map { meta, bed, bim, fam -> [meta.id, meta, bed, bim, fam] }
        .join(
            ch_headerless_phenotypes.filter { meta, _phenotype, _quant_covariates, _cat_covariates -> selects_direct.call(meta) }.map { meta, phenotype, quant_covariates, cat_covariates -> [meta.id, phenotype, quant_covariates, cat_covariates] },
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_ldak_weights.filter { meta, _weights -> selects_direct.call(meta) }.map { meta, weights -> [meta.id, weights] },
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .flatMap { _analysis_id, meta, bed, bim, fam, phenotype, quant_covariates, cat_covariates, weights ->
            meta.heritability_methods
                .findAll { method -> method in direct_methods }
                .collect { method ->
                    [meta, native_modes[getMethodCapability(method, 'estimator_family')], bed, bim, fam, phenotype, quant_covariates, cat_covariates, weights ?: []]
                }
        }
        .branch { _meta, mode, _bed, _bim, _fam, _phenotype, _quant_covariates, _cat_covariates, _weights ->
            he: mode == 'he'
            pcgc: mode == 'pcgc'
        }

    // The genotype bundle is cohort-scoped, so its role identity is the cohort. That is what makes the atom's
    // two-identity tag name the analysis and the bundle rather than repeating the analysis twice. The bundle
    // meta is tag-only: every output keys on the focal analysis meta.
    def ch_he_inputs = ch_direct_inputs.he.multiMap { meta, _mode, bed, bim, fam, phenotype, quant_covariates, cat_covariates, weights ->
        pheno: [meta, phenotype]
        genotypes: [[id: meta.cohort], bed, bim, fam, meta.method_options.ldak.power]
        weights: [meta, weights]
        qcovar: [meta, quant_covariates]
        covar: [meta, cat_covariates]
    }

    LDAK_FASTHE(
        ch_he_inputs.pheno,
        ch_he_inputs.genotypes,
        ch_he_inputs.weights,
        ch_he_inputs.qcovar,
        ch_he_inputs.covar,
    )

    def ch_pcgc_inputs = ch_direct_inputs.pcgc.multiMap { meta, _mode, bed, bim, fam, phenotype, quant_covariates, cat_covariates, weights ->
        pheno: [meta, phenotype, meta.population_prevalence]
        genotypes: [[id: meta.cohort], bed, bim, fam, meta.method_options.ldak.power]
        weights: [meta, weights]
        qcovar: [meta, quant_covariates]
        covar: [meta, cat_covariates]
    }

    LDAK_FASTPCGC(
        ch_pcgc_inputs.pheno,
        ch_pcgc_inputs.genotypes,
        ch_pcgc_inputs.weights,
        ch_pcgc_inputs.qcovar,
        ch_pcgc_inputs.covar,
    )

    emit:
    ldak_fast_he_results    = LDAK_FASTHE.out.fasthe_results // channel: [ val(meta), path(fasthe) ], one per analysis selecting ldak_fast_he
    ldak_fast_pcgc_results  = LDAK_FASTPCGC.out.fastpcgc_results // channel: [ val(meta), path(fastpcgc) ], one per analysis selecting ldak_fast_pcgc
    ldak_fast_pcgc_marginal = LDAK_FASTPCGC.out.fastpcgc_marginal // channel: [ val(meta), path(fastpcgc_marginal) ], where LDAK wrote it
}
