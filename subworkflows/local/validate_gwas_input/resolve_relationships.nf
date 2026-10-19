include { normaliseCellValue ; tokenizeMethodSelector } from './manifest_contracts'
include { getMethodCapabilities ; getRelationshipMethodTokens } from './method_registry'
include { getMethodOptionDefaults       } from './method_options'
include { getGctaBivariatePrevalence ; resolvePairLdmsMatrixSettings } from './request_contracts'

def resolveRelationships(relationship_rows, relationship_columns, relationship_manifest, analyses_by_id, summaries_by_id, errors) {
    def line_by_relationship_id = [:]
    def gcta_primary_requests = []
    def summary_pair_primary_requests = []
    def seen_bindings = [:]
    def referenced_analyses = [] as Set
    def referenced_summaries = [] as Set
    relationship_rows.eachWithIndex { row, index ->
        def line = index + 2
        def relationship_meta = row[0]
        def cells = [relationship_columns, row[1..-1]].transpose().collectEntries()
        def relationship_id = relationship_meta.id
        def left_analysis_id = normaliseCellValue(relationship_meta.left_analysis_id)?.toString()
        def right_analysis_id = normaliseCellValue(relationship_meta.right_analysis_id)?.toString()
        def left_summary_statistics_id = normaliseCellValue(relationship_meta.left_summary_statistics_id)?.toString()
        def right_summary_statistics_id = normaliseCellValue(relationship_meta.right_summary_statistics_id)?.toString()
        def reject = { field, message ->
            def named = field instanceof List ? field : [field]
            def label = named.size() > 1
                ? "fields ${named.collect { name -> "'${name}'" }.join(', ')}"
                : "field '${named.first()}'"
            errors << "  - ${relationship_manifest} row ${line} (relationship_id '${relationship_id}'), ${label}: ${message}"
        }

        if (line_by_relationship_id.containsKey(relationship_id)) {
            reject.call('relationship_id', "duplicate relationship_id '${relationship_id}', already declared on row ${line_by_relationship_id[relationship_id]}")
        }
        else {
            line_by_relationship_id[relationship_id] = line
        }
        def methods = tokenizeMethodSelector(relationship_meta.relationship_methods)
        def unknown = methods.findAll { method -> !getRelationshipMethodTokens().contains(method) }.unique()
        if (unknown) {
            reject.call('relationship_methods', "unknown method${unknown.size() > 1 ? 's' : ''} ${unknown.collect { method -> "'${method}'" }.join(', ')}, accepted values are ${getRelationshipMethodTokens().collect { method -> "'${method}'" }.join(', ')}")
        }
        def repeated = methods.countBy { method -> method }.findAll { _method, count -> count > 1 }.keySet()
        if (repeated) {
            reject.call('relationship_methods', "method${repeated.size() > 1 ? 's' : ''} ${repeated.collect { method -> "'${method}'" }.join(', ')} listed more than once")
        }
        if (!methods) {
            reject.call('relationship_methods', 'row selects no pairwise method')
        }
        def capabilities = getMethodCapabilities()
        def analysis_methods = methods.findAll { method -> capabilities[method]?.endpoint_domain == 'analysis' }
        def summary_methods = methods.findAll { method -> capabilities[method]?.endpoint_domain == 'summary_statistics' }

        if ((left_analysis_id && !right_analysis_id) || (!left_analysis_id && right_analysis_id)) {
            reject.call(['left_analysis_id', 'right_analysis_id'], 'analysis endpoint slots must be populated together')
        }
        if ((left_summary_statistics_id && !right_summary_statistics_id) || (!left_summary_statistics_id && right_summary_statistics_id)) {
            reject.call(['left_summary_statistics_id', 'right_summary_statistics_id'], 'summary-statistics endpoint slots must be populated together')
        }
        if (analysis_methods && (!left_analysis_id || !right_analysis_id)) {
            reject.call(['left_analysis_id', 'right_analysis_id'], "selected method${analysis_methods.size() > 1 ? 's' : ''} ${analysis_methods.join(', ')} require two analysis endpoints")
        }
        if (summary_methods && (!left_summary_statistics_id || !right_summary_statistics_id)) {
            reject.call(['left_summary_statistics_id', 'right_summary_statistics_id'], "selected method${summary_methods.size() > 1 ? 's' : ''} ${summary_methods.join(', ')} require two summary-statistics endpoints")
        }

        def left_analysis = left_analysis_id ? analyses_by_id[left_analysis_id] : null
        def right_analysis = right_analysis_id ? analyses_by_id[right_analysis_id] : null
        def left_summary = left_summary_statistics_id ? summaries_by_id[left_summary_statistics_id] : null
        def right_summary = right_summary_statistics_id ? summaries_by_id[right_summary_statistics_id] : null
        if (left_analysis_id && !left_analysis) {
            reject.call('left_analysis_id', "undefined analysis_id '${left_analysis_id}'")
        }
        if (right_analysis_id && !right_analysis) {
            reject.call('right_analysis_id', "undefined analysis_id '${right_analysis_id}'")
        }
        if (left_summary_statistics_id && !left_summary) {
            reject.call('left_summary_statistics_id', "undefined summary_statistics_id '${left_summary_statistics_id}'")
        }
        if (right_summary_statistics_id && !right_summary) {
            reject.call('right_summary_statistics_id', "undefined summary_statistics_id '${right_summary_statistics_id}'")
        }
        if (left_analysis_id && right_analysis_id && left_analysis_id == right_analysis_id) {
            reject.call(['left_analysis_id', 'right_analysis_id'], "the exact same analysis endpoint '${left_analysis_id}' cannot occupy both sides")
        }
        if (left_summary_statistics_id && right_summary_statistics_id && left_summary_statistics_id == right_summary_statistics_id) {
            reject.call(['left_summary_statistics_id', 'right_summary_statistics_id'], "the exact same summary-statistics endpoint '${left_summary_statistics_id}' cannot occupy both sides")
        }
        if (left_analysis && left_summary && left_summary.producer_analysis_id != left_analysis_id) {
            reject.call(['left_analysis_id', 'left_summary_statistics_id'], "same-side correspondence cannot be proven: summary '${left_summary_statistics_id}' was not produced by analysis '${left_analysis_id}'")
        }
        if (right_analysis && right_summary && right_summary.producer_analysis_id != right_analysis_id) {
            reject.call(['right_analysis_id', 'right_summary_statistics_id'], "same-side correspondence cannot be proven: summary '${right_summary_statistics_id}' was not produced by analysis '${right_analysis_id}'")
        }

        if (summary_methods && [left_summary, right_summary].any { summary -> summary?.source_kind == 'meta-analysis-derived' && summary.ancestry == 'MULTI' }) {
            reject.call(['left_summary_statistics_id', 'right_summary_statistics_id'], 'current LDSC and LDAK summary methods do not support a MULTI derived summary')
        }
        def left_meta = left_analysis ? left_analysis[0] : left_summary
        def right_meta = right_analysis ? right_analysis[0] : right_summary
        if (left_meta && right_meta && left_meta.trait.toString() == right_meta.trait.toString()) {
            reject.call(
                ['left_analysis_id', 'right_analysis_id', 'left_summary_statistics_id', 'right_summary_statistics_id'],
                "declared trait_id '${left_meta.trait}' is equal on both sides; self-pairs are invalid",
            )
        }
        if (analysis_methods && left_analysis && right_analysis && left_analysis[0].cohort.toString() != right_analysis[0].cohort.toString()) {
            reject.call(['left_analysis_id', 'right_analysis_id'], "individual-level relationship methods require one cohort, but '${left_analysis_id}' uses '${left_analysis[0].cohort}' and '${right_analysis_id}' uses '${right_analysis[0].cohort}'")
        }

        // A method whose declared trait support excludes binary endpoints must fail before execution rather
        // than return an observed-scale number nobody asked for. GCTA bivariate REML keeps the binary and
        // mixed pair domain because it is the only relationship estimator here with an explicit prevalence
        // and liability-scale contract.
        def binary_endpoint = (left_meta && left_meta.is_binary) || (right_meta && right_meta.is_binary)
        def quantitative_only_methods = binary_endpoint
            ? methods.findAll { method -> capabilities[method]?.trait_support && !capabilities[method].trait_support.binary }
            : []
        if (quantitative_only_methods) {
            reject.call(
                'relationship_methods',
                "method${quantitative_only_methods.size() > 1 ? 's' : ''} ${quantitative_only_methods.join(', ')} support two quantitative endpoints only; select 'gcta_bivariate_reml' or 'gcta_bivariate_reml_ldms' for a binary or mixed pair",
            )
        }

        // Declared pair covariates must reach the estimator or the request must fail. GCTA 1.94.1 accepts
        // `--qcovar`/`--covar` on an `--HEreg-bivar` command line and silently ignores them, so a
        // covariate-bearing HE request cannot be honoured and is refused here rather than answered wrongly.
        // Validate the relationship request: the programme accepts this combination but cannot honour it.
        // Retire when GCTA reads --covar/--qcovar under --HEreg-bivar.
        def covariate_incapable_methods = cells.pair_quant_covariates || cells.pair_cat_covariates
            ? methods.findAll { method -> capabilities[method]?.supports_covariates == false }
            : []
        if (covariate_incapable_methods) {
            reject.call(
                ['pair_quant_covariates', 'pair_cat_covariates'],
                "method${covariate_incapable_methods.size() > 1 ? 's' : ''} ${covariate_incapable_methods.join(', ')} have no native covariate parameter and would ignore the declared pair covariates; remove them or select a covariate-capable method",
            )
        }

        def left_side = "${left_analysis_id ?: ''}\u0001${left_summary_statistics_id ?: ''}".toString()
        def right_side = "${right_analysis_id ?: ''}\u0001${right_summary_statistics_id ?: ''}".toString()
        def binding_key = [left_side, right_side].sort().join('\u0000')
        if (seen_bindings.containsKey(binding_key)) {
            reject.call(
                ['left_analysis_id', 'right_analysis_id', 'left_summary_statistics_id', 'right_summary_statistics_id'],
                "unordered endpoint binding duplicates relationship_id '${seen_bindings[binding_key].id}' on row ${seen_bindings[binding_key].line}",
            )
        }
        else {
            seen_bindings[binding_key] = [id: relationship_id, line: line]
        }
        [left_analysis_id, right_analysis_id].findAll { endpoint -> endpoint }.each { endpoint -> referenced_analyses << endpoint }
        [left_summary_statistics_id, right_summary_statistics_id].findAll { endpoint -> endpoint }.each { endpoint -> referenced_summaries << endpoint }

        if (!unknown && !repeated && methods && left_meta && right_meta) {
            methods.each { method ->
                def request_id = "${method}--${relationship_id}".toString()
                def capability = capabilities[method]
                def common_meta = [
                    id: request_id,
                    request_id: request_id,
                    relationship_id: relationship_id,
                    method: method,
                    relationship_methods: methods,
                    left_analysis_id: left_analysis_id,
                    right_analysis_id: right_analysis_id,
                    left_summary_statistics_id: left_summary_statistics_id,
                    right_summary_statistics_id: right_summary_statistics_id,
                    left_trait_id: left_meta.trait,
                    right_trait_id: right_meta.trait,
                    left_trait_type: left_meta.trait_type,
                    right_trait_type: right_meta.trait_type,
                    left_is_binary: left_meta.is_binary,
                    right_is_binary: right_meta.is_binary,
                    left_genome_build: left_meta.build,
                    right_genome_build: right_meta.build,
                    left_ancestry: left_meta.ancestry,
                    right_ancestry: right_meta.ancestry,
                    left_population_prevalence: left_meta.population_prevalence,
                    right_population_prevalence: right_meta.population_prevalence,
                    left_sample_prevalence: left_meta.sample_prevalence,
                    right_sample_prevalence: right_meta.sample_prevalence,
                ]
                if (capability.endpoint_domain == 'analysis' && left_analysis && right_analysis) {
                    def left_analysis_meta = left_analysis[0]
                    def right_analysis_meta = right_analysis[0]
                    def matrix_settings = capability.component_model == 'ld_maf_stratified'
                        ? resolvePairLdmsMatrixSettings(
                            [:],
                            getMethodOptionDefaults().gcta.subMap(['ld_score_region_kb', 'ld_bins', 'ldms_maf_edges']),
                        ) { option, reason -> reject.call('relationship_methods', "default ${option}: ${reason}") }
                        : [:]
                    // A pair inherits the left endpoint's cohort identity, so it inherits the declared
                    // genotype view identity the same way and by the same rule: present only when declared,
                    // never as a `null` that would move every existing pair's task hash.
                    def declared_view_id = left_analysis_meta.genotype_view_id
                        ? [genotype_view_id: left_analysis_meta.genotype_view_id]
                        : [:]
                    def pair_meta = common_meta + [
                        cohort: left_analysis_meta.cohort,
                        build: left_analysis_meta.build,
                        ancestry: left_analysis_meta.ancestry,
                        genotype_format: left_analysis_meta.genotype_format,
                        matrix_kind: capability.matrix_kind,
                        matrix_settings: matrix_settings,
                        reml_bivar_prevalence: getGctaBivariatePrevalence(left_analysis_meta, right_analysis_meta),
                    ] + declared_view_id
                    def payload = [pair_meta, left_analysis[1], cells.pair_quant_covariates, cells.pair_cat_covariates]
                    gcta_primary_requests << [
                        request_id: request_id,
                        method: method,
                        reference_family: null,
                        meta: pair_meta,
                        payload: payload,
                    ]
                }
                if (capability.endpoint_domain == 'summary_statistics' && left_summary && right_summary) {
                    summary_pair_primary_requests << [
                        request_id: request_id,
                        method: method,
                        reference_family: capability.reference_family,
                        meta: common_meta,
                    ]
                }
            }
        }
    }
    return [
        analysis_primary_requests: gcta_primary_requests,
        summary_primary_requests: summary_pair_primary_requests,
        referenced_analyses: referenced_analyses,
        referenced_summaries: referenced_summaries,
        line_by_relationship_id: line_by_relationship_id,
    ]
}
