// Validate linked cohort and analysis manifests as one pipeline input contract and emit canonical analyses.
// This workflow performs no tasks and therefore emits no versions.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// PLUGIN
include { samplesheetToList                } from 'plugin/nf-schema'

// LOCAL FUNCTIONS
include { validateUniqueSamplesheetHeaders } from './manifest_contracts'
include { validateRelationalInput          } from './resolve_input'

workflow VALIDATE_GWAS_INPUT {
    take:
    cohort_manifest // channel: val(cohort_manifest)
    analysis_manifest // channel: val(analysis_manifest)
    summary_statistics_manifest // channel: val(summary_statistics_manifest)
    relationship_manifest // channel: val(relationship_manifest)
    reference_catalog // channel: val(reference_catalog)
    method_options // channel: val(method_options)

    main:
    if (analysis_manifest && !cohort_manifest) {
        error("[nf-core/gwas] ERROR: --analysis_manifest '${analysis_manifest ?: ''}' was supplied but --cohort_manifest is missing; linked manifests require both")
    }
    if (cohort_manifest && !analysis_manifest) {
        error("[nf-core/gwas] ERROR: --cohort_manifest '${cohort_manifest}' was supplied but --analysis_manifest is missing; linked manifests require both")
    }
    if (!analysis_manifest && !summary_statistics_manifest) {
        error("[nf-core/gwas] ERROR: supply either linked --cohort_manifest/--analysis_manifest inputs, --summary_statistics_manifest, or both")
    }

    def cohort_schema = "${projectDir}/assets/schema_cohort_manifest.json"
    def analysis_schema = "${projectDir}/assets/schema_analysis_manifest.json"
    def summary_statistics_schema = "${projectDir}/assets/schema_summary_statistics_manifest.json"
    def relationship_schema = "${projectDir}/assets/schema_relationship_manifest.json"
    if (cohort_manifest) {
        validateUniqueSamplesheetHeaders(cohort_manifest, 'Cohort manifest')
        validateUniqueSamplesheetHeaders(analysis_manifest, 'Analysis manifest')
    }
    if (summary_statistics_manifest) {
        validateUniqueSamplesheetHeaders(summary_statistics_manifest, 'Summary-statistics manifest')
    }
    if (relationship_manifest) {
        validateUniqueSamplesheetHeaders(relationship_manifest, 'Relationship manifest')
    }

    def validated = validateRelationalInput(
        cohort_manifest ? samplesheetToList(cohort_manifest, cohort_schema) : [],
        analysis_manifest ? samplesheetToList(analysis_manifest, analysis_schema) : [],
        summary_statistics_manifest ? samplesheetToList(summary_statistics_manifest, summary_statistics_schema) : [],
        relationship_manifest ? samplesheetToList(relationship_manifest, relationship_schema) : [],
        cohort_manifest,
        analysis_manifest,
        summary_statistics_manifest,
        relationship_manifest,
        cohort_schema,
        analysis_schema,
        summary_statistics_schema,
        relationship_schema,
        reference_catalog,
        method_options,
    )
    def ch_analyses = channel.fromList(validated.analyses)
    def ch_summary_statistics = channel.fromList(validated.summary_statistics)
    def ch_relationships = channel.fromList(validated.relationships)
    def ch_unary_requests = channel.fromList(validated.unary_requests)
    def ch_pair_requests = channel.fromList(validated.pair_requests)
    def ch_meta_requests = channel.fromList(validated.meta_requests)

    emit:
    analyses           = ch_analyses // channel: [ val(meta), [ path(genotype_file), ... ], path(phenotype), path(quant_covariates), path(cat_covariates), path(kvik_extract), path(ldak_weights) ]
    summary_statistics = ch_summary_statistics // channel: [ val(meta), path(source) ]
    relationships      = ch_relationships // channel: [ val(meta), [ path(genotype_file), ... ], path(pair_quant_covariates), path(pair_cat_covariates) ]
    unary_requests     = ch_unary_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ]
    meta_requests      = ch_meta_requests // channel: [ val(meta), val(source_summary_statistics_ids) ]
    pair_requests      = ch_pair_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ]
}
