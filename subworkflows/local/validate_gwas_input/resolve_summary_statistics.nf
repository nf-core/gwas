include { getInternalSummaryMetadata ; getSummaryStatisticsId } from './identity_helpers'
include { normaliseCellValue ; tokenizeMethodSelector } from './manifest_contracts'
include { getMethodCapabilities ; getSummaryUnaryMethodTokens } from './method_registry'

def resolveSummaryStatistics(summary_statistics_rows, summary_statistics_columns, summary_statistics_manifest, analyses_by_id, generated_summaries_by_id, errors) {
    def summaries_by_id = [:]
    def line_by_summary_statistics_id = [:]
    def validated_external_summaries = []
    def declared_summaries = []
    def meta_rows = []
    def summary_unary_primary_requests = []
    summary_statistics_rows.eachWithIndex { row, index ->
        def line = index + 2
        def summary_meta = row[0]
        def cells = [summary_statistics_columns, row[1..-1]].transpose().collectEntries()
        def summary_statistics_id = summary_meta.id
        def reject = { field, message ->
            def named = field instanceof List ? field : [field]
            def label = named.size() > 1
                ? "fields ${named.collect { name -> "'${name}'" }.join(', ')}"
                : "field '${named.first()}'"
            errors << "  - ${summary_statistics_manifest} row ${line} (summary_statistics_id '${summary_statistics_id}'), ${label}: ${message}"
        }

        if (line_by_summary_statistics_id.containsKey(summary_statistics_id)) {
            reject.call('summary_statistics_id', "duplicate summary_statistics_id '${summary_statistics_id}', already declared on row ${line_by_summary_statistics_id[summary_statistics_id]}")
        }
        else {
            line_by_summary_statistics_id[summary_statistics_id] = line
        }
        def source = normaliseCellValue(cells.source)
        def source_format = normaliseCellValue(summary_meta.source_format)?.toString()
        def producer_analysis_id = normaliseCellValue(summary_meta.producer_analysis_id)?.toString()
        def producer_association_method = normaliseCellValue(summary_meta.producer_association_method)?.toString()
        def has_external_origin = source || source_format
        def has_internal_origin = producer_analysis_id || producer_association_method
        def has_meta_origin = ['source_summary_statistics_ids', 'meta_analysis_models', 'min_studies']
            .any { field -> normaliseCellValue(summary_meta[field]) != null }
        if (has_meta_origin) {
            if (has_external_origin || has_internal_origin) {
                reject.call('source_summary_statistics_ids', 'meta-analysis origin is mutually exclusive with external and producer fields')
            }
            meta_rows << [meta: summary_meta, line: line]
            return
        }
        if (has_external_origin && has_internal_origin) {
            reject.call(
                ['source', 'source_format', 'producer_analysis_id', 'producer_association_method'],
                'external source fields and pipeline-generated producer fields are mutually exclusive',
            )
        }
        if (!has_external_origin && !has_internal_origin) {
            reject.call(
                ['source', 'producer_analysis_id'],
                'no summary-statistics origin is declared; supply a complete external source or producer analysis/method pair',
            )
        }
        if (has_external_origin && (!source || !source_format)) {
            reject.call(['source', 'source_format'], 'external origin requires both source fields')
        }
        if (has_internal_origin && (!producer_analysis_id || !producer_association_method)) {
            reject.call(['producer_analysis_id', 'producer_association_method'], 'pipeline-generated origin requires both producer fields')
        }
        def methods = tokenizeMethodSelector(summary_meta.heritability_methods)
        def unknown = methods.findAll { method -> !getSummaryUnaryMethodTokens().contains(method) }.unique()
        if (unknown) {
            reject.call('heritability_methods', "unknown unary summary method${unknown.size() > 1 ? 's' : ''} ${unknown.collect { method -> "'${method}'" }.join(', ')}, accepted values are ${getSummaryUnaryMethodTokens().collect { method -> "'${method}'" }.join(', ')}")
        }
        def repeated = methods.countBy { method -> method }.findAll { _method, count -> count > 1 }.keySet()
        if (repeated) {
            reject.call('heritability_methods', "method${repeated.size() > 1 ? 's' : ''} ${repeated.collect { method -> "'${method}'" }.join(', ')} listed more than once")
        }
        def resolved_meta = null
        if (has_internal_origin && producer_analysis_id && producer_association_method) {
            def producer = analyses_by_id[producer_analysis_id]
            if (!producer) {
                reject.call('producer_analysis_id', "undefined analysis_id '${producer_analysis_id}'")
            }
            else {
                if (!(producer_association_method in producer[0].association_methods)) {
                    reject.call(
                        'producer_association_method',
                        "analysis '${producer_analysis_id}' does not select association method '${producer_association_method}'",
                    )
                }
                def expected_id = getSummaryStatisticsId(producer_analysis_id, producer_association_method)
                if (summary_statistics_id != expected_id) {
                    reject.call(
                        'summary_statistics_id',
                        "pipeline-generated result must use deterministic identifier '${expected_id}'",
                    )
                }
                [
                    trait_id: normaliseCellValue(summary_meta.trait),
                    trait_type: normaliseCellValue(summary_meta.trait_type),
                    genome_build: normaliseCellValue(summary_meta.build),
                    ancestry: normaliseCellValue(summary_meta.ancestry),
                    source_method: normaliseCellValue(summary_meta.source_method),
                    source_release: normaliseCellValue(summary_meta.source_release),
                    population_prevalence: normaliseCellValue(summary_meta.population_prevalence),
                    sample_prevalence: normaliseCellValue(summary_meta.sample_prevalence),
                ].findAll { _field, value -> value != null }.each { field, _value ->
                    reject.call(field, "value is derived from producer analysis '${producer_analysis_id}'; leave this field blank")
                }
                resolved_meta = getInternalSummaryMetadata(producer[0], producer_association_method) + [
                    heritability_methods: methods,
                    access_constraints: normaliseCellValue(summary_meta.access_constraints),
                ]
            }
        }
        if (has_external_origin && source && source_format) {
            if (generated_summaries_by_id.containsKey(summary_statistics_id)) {
                reject.call(
                    'summary_statistics_id',
                    "external result collides with pipeline-generated result '${summary_statistics_id}'; choose a distinct external identity",
                )
            }
            def is_binary = summary_meta.trait_type == 'binary'
            def population_prevalence = normaliseCellValue(summary_meta.population_prevalence)
            def sample_prevalence = normaliseCellValue(summary_meta.sample_prevalence)
            if (!is_binary && population_prevalence != null) {
                reject.call('population_prevalence', "prevalence has no meaning on a quantitative trait")
            }
            if (!is_binary && sample_prevalence != null) {
                reject.call('sample_prevalence', "sample prevalence has no meaning on a quantitative trait")
            }
            resolved_meta = [
                id: summary_statistics_id,
                summary_statistics_id: summary_statistics_id,
                trait: summary_meta.trait,
                trait_id: summary_meta.trait,
                trait_type: summary_meta.trait_type,
                is_binary: is_binary,
                population_prevalence: population_prevalence,
                sample_prevalence: sample_prevalence,
                build: summary_meta.build,
                ancestry: summary_meta.ancestry,
                source_kind: 'external',
                source_format: source_format,
                source_method: summary_meta.source_method,
                source_release: normaliseCellValue(summary_meta.source_release),
                source_name: file(source).name,
                producer_analysis_id: null,
                producer_association_method: null,
                heritability_methods: methods,
                access_constraints: normaliseCellValue(summary_meta.access_constraints),
            ]
            validated_external_summaries << [resolved_meta, source]
        }
        if (resolved_meta) {
            if (!summaries_by_id.containsKey(summary_statistics_id)) {
                summaries_by_id[summary_statistics_id] = resolved_meta
                declared_summaries << resolved_meta
            }
            methods
                .findAll { method -> getSummaryUnaryMethodTokens().contains(method) }
                .each { method ->
                    def request_id = "${method}--${summary_statistics_id}".toString()
                    summary_unary_primary_requests << [
                        request_id: request_id,
                        method: method,
                        reference_family: getMethodCapabilities()[method].reference_family,
                        meta: resolved_meta + [
                            id: request_id,
                            request_id: request_id,
                            method: method,
                            summary_statistics_id: summary_statistics_id,
                        ],
                    ]
                }
        }
    }
    return [
        summaries: validated_external_summaries,
        declared_summaries: declared_summaries,
        meta_rows: meta_rows,
        summaries_by_id: summaries_by_id,
        unary_primary_requests: summary_unary_primary_requests,
        line_by_summary_statistics_id: line_by_summary_statistics_id,
    ]
}
