include { getMethodCapabilities               } from './method_registry'
include { getMethodCapability                 } from './method_registry'
include { getMethodTokensWithCapabilities     } from './method_registry'
include { getMethodOptionDefaults ; validateMphOptionValues } from './method_options'
include { validateAnalysisPairNativeArgumentTokens ; validateSummaryNativeArgumentTokens } from './native_option_policy'

def requestResourceTuple(meta, bundle) {
    return [
        meta,
        bundle && bundle.family == 'ldsc' ? bundle.resources.hapmap3_snplist : [],
        bundle && bundle.family == 'ldsc' ? bundle.resources.reference_ld_scores : [],
        bundle && bundle.family == 'ldsc' ? bundle.resources.regression_weights : [],
        bundle && bundle.family == 'ldak' ? bundle.resources.tagging_file : [],
    ]
}

def resolveRequestNamespace(method_options, document, namespace, primary_requests, reference_bundles) {
    def declared = method_options && document instanceof Map ? (document[namespace] ?: [:]) : [:]
    if (!(declared instanceof Map)) {
        error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}': expected an object")
    }
    def primaries = primary_requests.collectEntries { primary -> [(primary.request_id): primary] }
    def named = [:]
    declared.each { request_id, options ->
        if (!(options instanceof Map)) {
            error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request_id}': expected an option object")
        }
        if (!primaries.containsKey(request_id)) {
            def primary_request_id = options.primary_request_id
            def request_name = options.request_name
            if (!primary_request_id || !request_name) {
                error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request_id}': request is neither a selected deterministic primary nor a named addition with primary_request_id and request_name")
            }
            if (!primaries.containsKey(primary_request_id)) {
                error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request_id}': primary_request_id '${primary_request_id}' is not selected")
            }
            if (!(request_name instanceof String) || !(request_name ==~ /^[A-Za-z0-9_.+-]+$/)) {
                error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request_id}', option 'request_name': expected one non-empty identifier token")
            }
            def expected_id = "${primary_request_id}--${request_name}".toString()
            if (request_id != expected_id) {
                error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request_id}': named addition must use deterministic identifier '${expected_id}'")
            }
            named[request_id] = primaries[primary_request_id] + [
                request_id: request_id,
                primary_request_id: primary_request_id,
                request_name: request_name,
                is_primary: false,
            ]
        }
        else if (options.containsKey('primary_request_id') || options.containsKey('request_name')) {
            error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request_id}': deterministic primary requests cannot declare primary_request_id or request_name")
        }
    }

    def resolved = []
    (primary_requests.collect { primary -> primary + [primary_request_id: primary.request_id, request_name: null, is_primary: true] } + named.values()).each { request ->
        def options = declared[request.request_id] ?: [:]
        def accepted = request.reference_family
            ? ['reference_bundle_id', 'native_args', 'primary_request_id', 'request_name']
            : ['native_args', 'primary_request_id', 'request_name']
        // The LD- and MAF-stratified plan settings are accepted by the component model the estimator fits,
        // not by the matrix kind that happens to carry it today. That is what lets a second tool's stratified
        // estimator share the plan without this list learning its matrix kind.
        def fits_ld_maf_strata = getMethodCapability(request.method, 'component_model') == 'ld_maf_stratified'
        if (fits_ld_maf_strata) {
            accepted += ['ld_score_region_kb', 'ld_bins', 'ldms_maf_edges']
        }
        // The curated MPH controls are admitted by the estimator's option family, not by the namespace, so a
        // GCTA pair request in the same document still rejects them as unknown. `option_family` is optional in
        // the capability contract, so it is read with the optional accessor.
        def option_family = getMethodCapabilities()[request.method]?.option_family
        if (option_family == 'mph') {
            accepted += ['iterations', 'tolerance', 'random_vectors', 'seed', 'save_memory']
        }
        def unknown = options.keySet().findAll { option -> !(option in accepted) }
        if (unknown) {
            error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request.request_id}', option '${unknown.first()}': unknown option; accepted options are ${accepted.join(', ')}")
        }
        def request_options = [:]
        if (option_family == 'mph') {
            request_options = [mph: validateMphOptionValues(options, getMethodOptionDefaults().mph) { option, reason ->
                error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request.request_id}', option '${option}': ${reason}")
            }]
        }
        def native_args = request.reference_family
            ? validateSummaryNativeArgumentTokens(method_options, namespace, request.request_id, request.method, options.native_args ?: [])
            : validateAnalysisPairNativeArgumentTokens(method_options, request.request_id, request.method, options.native_args ?: [])
        def matrix_settings = request.meta.matrix_settings ?: [:]
        if (fits_ld_maf_strata) {
            def fail = { option, reason ->
                error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request.request_id}', option '${option}': ${reason}")
            }
            matrix_settings = resolvePairLdmsMatrixSettings(options, matrix_settings, fail)
        }
        def bundle = null
        if (request.reference_family) {
            def reference_bundle_id = options.reference_bundle_id
            if (!(reference_bundle_id instanceof String) || !reference_bundle_id.trim()) {
                error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request.request_id}', option 'reference_bundle_id': every summary-statistics request must explicitly select a reference bundle")
            }
            bundle = reference_bundles[reference_bundle_id]
            if (!bundle) {
                error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request.request_id}', option 'reference_bundle_id': undefined reference bundle '${reference_bundle_id}'")
            }
            if (bundle.family != request.reference_family) {
                error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request.request_id}', option 'reference_bundle_id': method '${request.method}' requires family '${request.reference_family}', but '${reference_bundle_id}' belongs to '${bundle.family}'")
            }
            // A one-category tagging model is what makes the SumCors component model single-component; the
            // registry names the pairwise LDAK summary estimator rather than this seam repeating its token.
            def ldak_summary_pairs = getMethodTokensWithCapabilities([domain: 'pairwise', reference_family: 'ldak'])
            if (request.method in ldak_summary_pairs && !(bundle.model in ['LDAK-Thin', 'Uniform-GCTA', 'Human-Default'])) {
                error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request.request_id}', option 'reference_bundle_id': first-release LDAK SumCors supports models 'LDAK-Thin', 'Uniform-GCTA' and 'Human-Default', but bundle '${reference_bundle_id}' declares '${bundle.model}'")
            }
        }
        def resolved_meta = request.meta + [
            id: request.request_id,
            request_id: request.request_id,
            primary_request_id: request.primary_request_id,
            request_name: request.request_name,
            is_primary: request.is_primary,
            native_args: native_args,
            request_options: request_options,
            matrix_settings: matrix_settings,
            reference_bundle_id: bundle ? bundle.id : null,
            reference_family: bundle ? bundle.family : null,
            reference_metadata: bundle
                ? bundle.subMap(['id', 'family', 'genome_build', 'ancestry', 'variant_id_system', 'model', 'role_names'])
                : null,
            reference_validation: bundle ? 'structural_availability_only' : null,
        ]
        resolved << [request: request, meta: resolved_meta, bundle: bundle]
    }
    return resolved
}

def resolvePairLdmsMatrixSettings(options, defaults, fail) {
    def ld_score_region_kb = options.containsKey('ld_score_region_kb') ? options.ld_score_region_kb : defaults.ld_score_region_kb
    def ld_bins = options.containsKey('ld_bins') ? options.ld_bins : defaults.ld_bins
    def ldms_maf_edges = options.containsKey('ldms_maf_edges')
        ? options.ldms_maf_edges
        : defaults.containsKey('maf_edges')
            ? defaults.maf_edges
            : defaults.ldms_maf_edges

    if (!(ld_score_region_kb instanceof Number) || ld_score_region_kb < 1 || ld_score_region_kb != ld_score_region_kb.toInteger()) {
        fail.call('ld_score_region_kb', 'expected one positive integer')
    }
    if (!(ld_bins instanceof Number) || ld_bins < 1 || ld_bins != ld_bins.toInteger()) {
        fail.call('ld_bins', 'expected one positive integer')
    }
    if (!(ldms_maf_edges instanceof List) || ldms_maf_edges.size() < 2 || !ldms_maf_edges.every { edge -> edge instanceof Number && edge >= 0 && edge <= 0.5 }) {
        fail.call('ldms_maf_edges', 'expected a numeric boundary list spanning 0 to 0.5')
    }
    if (ldms_maf_edges.first() != 0 || ldms_maf_edges.last() != 0.5) {
        fail.call('ldms_maf_edges', 'boundaries must start at 0 and end at 0.5')
    }
    if ((1..<ldms_maf_edges.size()).any { index -> ldms_maf_edges[index] <= ldms_maf_edges[index - 1] }) {
        fail.call('ldms_maf_edges', 'boundaries must be strictly increasing')
    }
    return [
        ld_score_region_kb: ld_score_region_kb.toInteger(),
        ld_bins: ld_bins.toInteger(),
        maf_edges: ldms_maf_edges,
    ]
}

def getGctaBivariatePrevalence(left_meta, right_meta) {
    if (left_meta.is_binary && right_meta.is_binary) {
        return left_meta.population_prevalence != null && right_meta.population_prevalence != null
            ? [left_meta.population_prevalence, right_meta.population_prevalence]
            : []
    }
    if (left_meta.is_binary && left_meta.population_prevalence != null) {
        return [left_meta.population_prevalence]
    }
    if (right_meta.is_binary && right_meta.population_prevalence != null) {
        return [right_meta.population_prevalence]
    }
    return []
}
