include { normaliseCellValue ; tokenizeMethodSelector } from './manifest_contracts'
include { getSummarySetMethodTokens ; getSummaryUnaryMethodTokens ; getMethodCapabilities } from './method_registry'

// Resolve the manifest's N-way producer only after every base summary has entered the registry.
// Programme counts, numerical fit availability and native statistical output are not ingress constraints.
def resolveMetaAnalysisRequests(meta_rows, summaries, generated_summaries_by_id, options_document, method_options, manifest, errors) {
    def requests = []
    def referenced_summaries = [] as Set
    def declared_ids = meta_rows.collect { row -> row.meta.id } as Set
    def options_by_id = options_document.meta_requests ?: [:]
    options_by_id.keySet().findAll { id -> !(id in declared_ids) }.each { id ->
        errors << "  - Method-options document '${method_options}', namespace 'meta_requests', request_id '${id}': no meta-analysis row declares this result"
    }
    meta_rows.each { row ->
        def meta = row.meta
        def id = meta.id
        def errors_before = errors.size()
        def reject = { field, reason ->
            errors << "  - ${manifest} row ${row.line} (summary_statistics_id '${id}'), field '${field}': ${reason}"
        }
        def source_ids = tokenizeMethodSelector(meta.source_summary_statistics_ids)
        def models = tokenizeMethodSelector(meta.meta_analysis_models) ?: ['fixed']
        if (summaries.summaries_by_id.containsKey(id) || generated_summaries_by_id.containsKey(id)) {
            reject.call('summary_statistics_id', 'meta-analysis identity collides with a base summary identity')
        }
        if (source_ids.size() < 2) {
            reject.call('source_summary_statistics_ids', 'declare at least two base summary identifiers')
        }
        if (source_ids.unique(false).size() != source_ids.size()) {
            reject.call('source_summary_statistics_ids', 'a parent summary is listed more than once')
        }
        source_ids.each { source_id ->
            if (source_id == id) {
                reject.call('source_summary_statistics_ids', 'a result cannot be its own parent')
            }
            else if (source_id in declared_ids) {
                reject.call('source_summary_statistics_ids', "parent '${source_id}' is meta-analysis-derived; parents must be declared base summaries")
            }
            else if (!summaries.summaries_by_id.containsKey(source_id)) {
                reject.call('source_summary_statistics_ids', "undefined declared base summary '${source_id}'")
            }
        }
        def unknown_models = models.findAll { model -> !(model in getSummarySetMethodTokens()) }
        if (unknown_models) {
            reject.call('meta_analysis_models', "unknown model token(s): ${unknown_models.join(', ')}")
        }
        if (models.unique(false).size() != models.size()) {
            reject.call('meta_analysis_models', 'a model is listed more than once')
        }
        if (!models.contains('fixed')) {
            reject.call('meta_analysis_models', 'fixed is required alongside every optional model')
        }
        def min_studies = normaliseCellValue(meta.min_studies)
        min_studies = min_studies == null ? 2 : min_studies
        // This is the declared N-way result selection, independently of every programme's native limits.
        if (!(min_studies instanceof Number) || min_studies != min_studies.toInteger() || min_studies < 2 || min_studies > source_ids.size()) {
            reject.call('min_studies', 'expected an integer from 2 through the declared parent count')
        }
        def options = options_by_id.containsKey(id) ? options_by_id[id] : [:]
        def axes = null
        if (!(options instanceof Map)) {
            reject.call('method_options.meta_requests', 'expected an option object keyed by the result identifier')
        }
        else {
            def unknown_options = options.keySet().findAll { option -> option != 'mrmega' }
            if (unknown_options) {
                reject.call('method_options.meta_requests', "unknown options: ${unknown_options.join(', ')}")
            }
            if (options.containsKey('mrmega') && !models.contains('mrmega')) {
                reject.call('method_options.meta_requests.mrmega', 'mrmega must be selected on this summary row')
            }
            if (models.contains('mrmega')) {
                def mrmega = options.mrmega
                if (!(mrmega instanceof Map) || mrmega.keySet() != ['axes'].toSet()) {
                    reject.call('method_options.meta_requests.mrmega', 'declare an object containing axes')
                }
                else if (!(mrmega.axes instanceof Number) || mrmega.axes != mrmega.axes.toInteger()) {
                    reject.call('method_options.meta_requests.mrmega.axes', 'expected an integer')
                }
                else {
                    axes = mrmega.axes.toInteger()
                }
            }
        }
        ['trait', 'trait_type', 'build', 'ancestry', 'source_method'].each { field ->
            if (normaliseCellValue(meta[field]) != null) {
                reject.call(field, 'value is derived from the declared parent summaries; leave this field blank')
            }
        }
        def methods = tokenizeMethodSelector(meta.heritability_methods)
        def unknown_methods = methods.findAll { method -> !(method in getSummaryUnaryMethodTokens()) }
        if (unknown_methods) {
            reject.call('heritability_methods', "unknown unary summary methods: ${unknown_methods.join(', ')}")
        }
        if (methods.unique(false).size() != methods.size()) {
            reject.call('heritability_methods', 'a method is listed more than once')
        }
        if (errors.size() != errors_before) {
            return
        }
        def parents = source_ids.collect { source_id -> summaries.summaries_by_id[source_id] }
        ['trait', 'trait_type', 'build'].each { field ->
            if (parents.collect { parent -> parent[field] }.unique().size() != 1) {
                reject.call('source_summary_statistics_ids', "parents declare incompatible ${field} values")
            }
        }
        def ancestries = parents.collect { parent -> parent.ancestry }.unique().sort()
        if (ancestries.size() > 1 && !models.contains('mrmega')) {
            reject.call('meta_analysis_models', 'parents with different ancestry labels require explicit mrmega selection')
        }
        def ancestry = ancestries.size() == 1 ? ancestries.first() : 'MULTI'
        if (ancestry == 'MULTI' && methods) {
            reject.call('heritability_methods', 'current LDSC and LDAK summary methods do not support a MULTI derived summary')
        }
        def population_prevalence = normaliseCellValue(meta.population_prevalence)
        def sample_prevalence = normaliseCellValue(meta.sample_prevalence)
        def inherited_prevalences = parents.collect { parent -> parent.population_prevalence }.findAll { value -> value != null }.unique()
        if (population_prevalence == null && inherited_prevalences.size() == 1) {
            population_prevalence = inherited_prevalences.first()
        }
        if (parents.first().trait_type != 'binary' && (population_prevalence != null || sample_prevalence != null)) {
            reject.call('population_prevalence, sample_prevalence', 'prevalence has no meaning on a quantitative trait')
        }
        if (errors.size() != errors_before) {
            return
        }
        def result_meta = [
            id: id,
            request_id: id,
            summary_statistics_id: id,
            method: 'common_variant_meta_analysis',
            source_kind: 'meta-analysis-derived',
            source_format: 'gwaslab',
            source_method: 'fixed',
            source_release: normaliseCellValue(meta.source_release),
            source_name: "${id}.gwaslab.tsv.gz".toString(),
            producer_analysis_id: null,
            producer_association_method: null,
            source_summary_statistics_ids: source_ids,
            source_ancestries: parents.collect { parent -> parent.ancestry },
            ancestry_components: ancestries,
            meta_analysis_models: models,
            min_studies: min_studies.toInteger(),
            axes: axes,
            trait: parents.first().trait,
            trait_id: parents.first().trait,
            trait_type: parents.first().trait_type,
            is_binary: parents.first().is_binary,
            build: parents.first().build,
            ancestry: ancestry,
            population_prevalence: population_prevalence,
            sample_prevalence: sample_prevalence,
            heritability_methods: methods,
            access_constraints: (parents.collect { parent -> parent.access_constraints } + normaliseCellValue(meta.access_constraints))
                .findAll { value -> value }.unique().join('; ') ?: null,
        ]
        summaries.summaries_by_id[id] = result_meta
        summaries.declared_summaries << result_meta
        methods.each { method ->
            def request_id = "${method}--${id}".toString()
            summaries.unary_primary_requests << [
                request_id: request_id,
                method: method,
                reference_family: getMethodCapabilities()[method].reference_family,
                meta: result_meta + [id: request_id, request_id: request_id, method: method],
            ]
        }
        requests << [result_meta, source_ids]
        referenced_summaries.addAll(source_ids)
    }
    return [requests: requests, referenced_summaries: referenced_summaries]
}
