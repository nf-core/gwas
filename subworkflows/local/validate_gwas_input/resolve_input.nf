include { getSamplesheetPositionalColumns } from './manifest_contracts'
include {
    readMethodOptionsDocument ;
    validateMethodOptions
} from './method_options'
include {
    getPrevalenceConsumers ;
    resolveAnalyses
} from './resolve_analyses'
include { resolveCohorts                  } from './resolve_cohorts'
include { readReferenceCatalog            } from './resolve_references'
include { resolveRelationships            } from './resolve_relationships'
include { resolveMetaAnalysisRequests } from './resolve_meta_analysis'
include { resolveSummaryStatistics        } from './resolve_summary_statistics'
include {
    requestResourceTuple ;
    resolveRequestNamespace
} from './request_contracts'

// Validate the lists as a linked contract and construct the canonical tuple consumed by GWAS.
def validateRelationalInput(cohort_rows, analysis_rows, summary_statistics_rows, relationship_rows, cohort_manifest, analysis_manifest, summary_statistics_manifest, relationship_manifest, cohort_schema, analysis_schema, summary_statistics_schema, relationship_schema, reference_catalog, method_options = null) {
    def errors = []
    def options_document = readMethodOptionsDocument(method_options)
    def method_options_by_analysis = validateMethodOptions(method_options, analysis_rows, options_document)
    def reference_bundles = readReferenceCatalog(reference_catalog)
    def prevalence_consumers = getPrevalenceConsumers(relationship_rows, summary_statistics_rows)

    def cohorts = resolveCohorts(
        cohort_rows,
        getSamplesheetPositionalColumns(cohort_schema),
        cohort_manifest,
        errors,
    )
    def analyses = resolveAnalyses(
        analysis_rows,
        getSamplesheetPositionalColumns(analysis_schema),
        analysis_manifest,
        cohort_manifest,
        cohorts.cohorts_by_id,
        relationship_rows,
        prevalence_consumers,
        method_options_by_analysis,
        errors,
    )
    def summaries = resolveSummaryStatistics(
        summary_statistics_rows,
        getSamplesheetPositionalColumns(summary_statistics_schema),
        summary_statistics_manifest,
        analyses.analyses_by_id,
        analyses.generated_summaries_by_id,
        errors,
    )
    def meta_analysis = resolveMetaAnalysisRequests(
        summaries.meta_rows, summaries, analyses.generated_summaries_by_id,
        options_document, method_options, summary_statistics_manifest, errors,
    )
    def relationships = resolveRelationships(
        relationship_rows,
        getSamplesheetPositionalColumns(relationship_schema),
        relationship_manifest,
        analyses.analyses_by_id,
        summaries.summaries_by_id,
        errors,
    )

    analyses.analyses.each { row ->
        def meta = row[0]
        if (!meta.association_methods && !meta.heritability_methods && !(meta.id in relationships.referenced_analyses)) {
            def line = analyses.line_by_analysis_id[meta.id]
            errors << "  - ${analysis_manifest} row ${line} (analysis_id '${meta.id}'), fields 'association_methods', 'heritability_methods': row selects no unary method and is not referenced by any relationship; populate a method, reference it from a relationship or remove the row"
        }
    }
    summaries.declared_summaries.each { meta ->
        if (meta.source_kind != 'meta-analysis-derived' && !meta.heritability_methods && !(meta.summary_statistics_id in relationships.referenced_summaries) && !(meta.summary_statistics_id in meta_analysis.referenced_summaries)) {
            def line = summaries.line_by_summary_statistics_id[meta.summary_statistics_id]
            errors << "  - ${summary_statistics_manifest} row ${line} (summary_statistics_id '${meta.summary_statistics_id}'), field 'heritability_methods': result selects no unary method and is not referenced by any relationship"
        }
    }
    if (errors) {
        error("[nf-core/gwas] ERROR: Validation of linked manifests failed!\n\n${errors.join('\n')}\n")
    }

    def resolved_unary = resolveRequestNamespace(
        method_options,
        options_document,
        'unary_requests',
        summaries.unary_primary_requests,
        reference_bundles,
    )
    def resolved_pair = resolveRequestNamespace(
        method_options,
        options_document,
        'pair_requests',
        relationships.analysis_primary_requests + relationships.summary_primary_requests,
        reference_bundles,
    )
    def resolved_relationships = resolved_pair
        .findAll { resolved -> !resolved.request.reference_family }
        .collect { resolved ->
            def payload = resolved.request.payload
            [resolved.meta, payload[1], payload[2], payload[3]]
        }
    def unary_requests = resolved_unary.collect { resolved -> requestResourceTuple(resolved.meta, resolved.bundle) }
    def pair_requests = resolved_pair
        .findAll { resolved -> resolved.request.reference_family }
        .collect { resolved -> requestResourceTuple(resolved.meta, resolved.bundle) }

    return [
        analyses: analyses.analyses,
        summary_statistics: summaries.summaries,
        relationships: resolved_relationships,
        unary_requests: unary_requests,
        pair_requests: pair_requests,
        meta_requests: meta_analysis.requests,
    ]
}
