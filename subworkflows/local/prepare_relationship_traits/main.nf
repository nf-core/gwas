// Prepare each declared relationship's ordered two-trait table and its relationship-owned covariates exactly
// once, whichever individual-level pair estimators select it. A relationship is the scientific pair: the
// ordered full-union phenotype table and the normalised pair covariates are properties of the two endpoints,
// not of the estimator that consumes them, so every method adapter serialises this one artifact rather than
// resolving the endpoints a second time and defining a second endpoint sample set.
//
// This owns relationship de-duplication and endpoint resolution against the prepared unary phenotypes, which
// is relational-input policy, so it is pipeline-owned rather than a component-library submission candidate.
// It reads no params, no workflow and no projectDir, and every constituent process reports directly to the
// run-wide versions topic, so it emits no versions.

// MODULE: Local to the pipeline
include { PREPARE_BIVARIATE_TRAITS } from '../../../modules/local/prepare_bivariate_traits/main'

workflow PREPARE_RELATIONSHIP_TRAITS {
    take:
    ch_relationships // channel: [ val(meta), [ path(genotype_file), ... ], path(pair_quant_covariates), path(pair_cat_covariates) ], one validated row per selected method per relationship
    ch_prepared_phenotypes // channel: [ val(meta), path(phenotype) ], the headerless phenotype of every analysis unit, keyed one-to-one on the analysis meta

    main:

    // The prepared pair is reusable work shared by every estimator that selects the relationship, so its
    // identity is built from an explicit list of relationship-defining keys rather than from whichever
    // request's row happened to survive `unique`. `method`, `native_args`, `matrix_kind`, `matrix_settings`
    // and `request_options` are deliberately absent: they are request state, and this map is the `val(meta)`
    // input of PREPARE_BIVARIATE_TRAITS and therefore part of its task hash. Carrying them would let a change
    // of one estimator's seed invalidate a table that estimator does not define.
    def ch_relationship_definitions = ch_relationships
        .map { meta, genotype_files, pair_quant_covariates, pair_cat_covariates ->
            [meta.relationship_id, buildRelationshipMeta(meta), genotype_files, pair_quant_covariates, pair_cat_covariates]
        }
        .unique { relationship_id, _meta, _genotype_files, _pair_quant_covariates, _pair_cat_covariates -> relationship_id }

    // `combine` is deliberate at the endpoint seams: one analysis unit may be an endpoint of several
    // relationships, so neither side is one-to-one until the two are joined back on the relationship.
    def ch_left_pair_phenotypes = ch_relationship_definitions
        .map { relationship_id, meta, _genotype_files, pair_quant_covariates, pair_cat_covariates ->
            [meta.left_analysis_id, relationship_id, meta, pair_quant_covariates ?: [], pair_cat_covariates ?: []]
        }
        .combine(
            ch_prepared_phenotypes.map { meta, phenotype -> [meta.id, phenotype] },
            by: 0
        )
        .map { _analysis_id, relationship_id, meta, pair_quant_covariates, pair_cat_covariates, phenotype ->
            [relationship_id, meta, phenotype, pair_quant_covariates, pair_cat_covariates]
        }

    def ch_right_pair_phenotypes = ch_relationship_definitions
        .map { relationship_id, meta, _genotype_files, _pair_quant_covariates, _pair_cat_covariates -> [meta.right_analysis_id, relationship_id] }
        .combine(
            ch_prepared_phenotypes.map { meta, phenotype -> [meta.id, phenotype] },
            by: 0
        )
        .map { _analysis_id, relationship_id, phenotype -> [relationship_id, phenotype] }

    def ch_pair_trait_inputs = ch_left_pair_phenotypes
        .join(ch_right_pair_phenotypes, failOnDuplicate: true, failOnMismatch: true)
        .map { _relationship_id, meta, left_phenotype, pair_quant_covariates, pair_cat_covariates, right_phenotype ->
            [meta, left_phenotype, right_phenotype, pair_quant_covariates, pair_cat_covariates]
        }

    //
    // MODULE: one ordered two-trait table and one normalised pair covariate set per relationship
    //
    PREPARE_BIVARIATE_TRAITS(ch_pair_trait_inputs)

    // Both emits carry exactly one record per relationship. The optional covariate outputs are folded in with
    // `remainder: true` and rendered as `[]`, so a covariate-free relationship still emits a record and the
    // consumer's join stays a guarded one-to-one rather than silently losing the relationship.
    def ch_pairs = PREPARE_BIVARIATE_TRAITS.out.phenotype
        .join(PREPARE_BIVARIATE_TRAITS.out.quant_covariates, remainder: true)
        .join(PREPARE_BIVARIATE_TRAITS.out.cat_covariates, remainder: true)
        .map { meta, phenotype, quant_covariates, cat_covariates ->
            [meta, phenotype, quant_covariates ?: [], cat_covariates ?: []]
        }

    def ch_named_covariates = PREPARE_BIVARIATE_TRAITS.out.phenotype
        .join(PREPARE_BIVARIATE_TRAITS.out.named_quant_covariates, remainder: true)
        .join(PREPARE_BIVARIATE_TRAITS.out.named_cat_covariates, remainder: true)
        .map { meta, _phenotype, quant_covariates, cat_covariates ->
            [meta, quant_covariates ?: [], cat_covariates ?: []]
        }

    emit:
    pairs = ch_pairs // channel: [ val(relationship_meta), path(phenotype), path(quant_covariates), path(cat_covariates) ], headerless covariates, [] when absent
    named_covariates = ch_named_covariates // channel: [ val(relationship_meta), path(quant_covariates), path(cat_covariates) ], headered covariates, [] when absent
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// The keys that define the scientific pair, and nothing else. Two request rows of one relationship that differ
// only in their estimator state must produce the identical map, because that map is the preparation task's
// identity and the preparation is shared.
def getRelationshipDefinitionKeys() {
    return [
        'relationship_id',
        'cohort',
        'build',
        'ancestry',
        'genotype_format',
        'left_analysis_id',
        'right_analysis_id',
        'left_trait_id',
        'right_trait_id',
        'left_trait_type',
        'right_trait_type',
        'left_is_binary',
        'right_is_binary',
        'left_genome_build',
        'right_genome_build',
        'left_ancestry',
        'right_ancestry',
    ]
}

def buildRelationshipMeta(meta) {
    return meta.subMap(getRelationshipDefinitionKeys()) + [
        id: meta.relationship_id,
        request_id: meta.relationship_id,
    ]
}
