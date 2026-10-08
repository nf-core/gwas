include { getMethodCapabilities ; getMethodTokensWithCapabilities } from './method_registry'
include { tokenizeMethodSelector          } from './manifest_contracts'

// The accepted method families, their accepted option names and every resolved default are this one map.
// Family acceptance, per-family key validation and the resolved shape are all derived from it, so adding a
// family is adding an entry here rather than editing three parallel lists. Families are keyed alphabetically.
//
// A `null` default always means "pass nothing and let the native default stand", never "zero" or "off": the
// three `fast_*` LDAK controls omit `--repetitions`, `--num-blocks` and `--random-seed` respectively, which is
// what makes an unset option reproduce the tool's own documented behaviour rather than a pipeline opinion.
def getMethodOptionDefaults() {
    return [
        gcta: [
            grm_maf: null,
            grm_extract: [],
            reml_no_constrain: false,
            sparse_cutoff: 0.05,
            ld_score_region_kb: 200,
            ld_bins: 4,
            ldms_maf_edges: [0, 0.01, 0.05, 0.2, 0.5],
        ],
        ldak: [
            model: 'human_default',
            power: -0.25,
            weights_policy: 'equal',
            weights: [],
            relatedness_filter: false,
            kvik_step1_subset: 'all',
            predictor_extract: [],
            kvik_step2_keep: [],
            fast_repetitions: null,
            fast_num_blocks: null,
            fast_seed: null,
        ],
        mph: [
            iterations: null,
            tolerance: null,
            random_vectors: null,
            seed: null,
            save_memory: false,
        ],
        regenie: [
            step1_bsize: 1000,
            firth: true,
            firth_approx: true,
            firth_p_threshold: 0.01,
            min_mac: null,
        ],
    ]
}

// Operational tuning that a researcher might reasonably try to place in the scientific option document. Naming
// it per family lets the diagnostic say "run/profile configuration" instead of "unknown option".
def getMethodOptionOperationalKeys() {
    return [
        gcta: ['gcta_grm_parts'],
        ldak: ['threads', 'jobs', 'partitions'],
        mph: ['threads', 'num_threads'],
        regenie: ['step2_bsize', 'step1_mode', 'step1_jobs', 'lowmem'],
    ]
}

def validateMethodOptionKeys(analysis_id, family, options, accepted, operational, fail) {
    def unknown = options.keySet().findAll { option -> !(option in accepted) }
    if (unknown) {
        def option = unknown.first()
        def reason = option in operational
            ? 'operational tuning must be supplied through run/profile configuration'
            : "unknown option; accepted ${family.toUpperCase()} options are ${accepted.join(', ')}"
        if (family == 'gcta' && option == 'gcta_grm_parts') {
            reason = 'partition count is operational and must be supplied through run/profile configuration'
        }
        fail.call(analysis_id, "${family}.${option}", reason)
    }
}

def resolveMethodResource(analysis_id, family, option, options, fail) {
    if (!options.containsKey(option)) {
        return []
    }
    def declared = options[option]
    if (!(declared instanceof String) || !declared.trim()) {
        fail.call(analysis_id, "${family}.${option}", 'expected a non-empty resource path string')
    }
    def resolved = file(declared)
    if (!resolved.exists()) {
        fail.call(analysis_id, "${family}.${option}", "resource path '${declared}' does not exist")
    }
    return resolved
}

def resolveRegenieMethodOptions(analysis_id, options, methods, defaults, fail) {
    def regenie_association = getMethodTokensWithCapabilities([domain: 'association', option_family: 'regenie'])
    if (options && !methods.association_methods.any { method -> method in regenie_association }) {
        fail.call(analysis_id, "regenie.${options.keySet().first()}", "analysis does not select 'regenie'")
    }

    def step1_bsize = options.containsKey('step1_bsize') ? options.step1_bsize : defaults.step1_bsize
    def firth = options.containsKey('firth') ? options.firth : defaults.firth
    def firth_approx = options.containsKey('firth_approx') ? options.firth_approx : defaults.firth_approx
    def firth_p_threshold = options.containsKey('firth_p_threshold') ? options.firth_p_threshold : defaults.firth_p_threshold
    def min_mac = options.containsKey('min_mac') ? options.min_mac : defaults.min_mac
    if (!(step1_bsize instanceof Number) || step1_bsize < 1 || step1_bsize != step1_bsize.toInteger()) {
        fail.call(analysis_id, 'regenie.step1_bsize', 'expected one positive integer')
    }
    if (!(firth instanceof Boolean)) {
        fail.call(analysis_id, 'regenie.firth', 'expected a boolean')
    }
    if (!(firth_approx instanceof Boolean)) {
        fail.call(analysis_id, 'regenie.firth_approx', 'expected a boolean')
    }
    if (options.containsKey('firth_approx') && firth_approx && !firth) {
        fail.call(analysis_id, 'regenie.firth_approx', "requires 'regenie.firth' to be true")
    }
    if (!(firth_p_threshold instanceof Number) || firth_p_threshold <= 0 || firth_p_threshold > 1) {
        fail.call(analysis_id, 'regenie.firth_p_threshold', 'expected a number greater than 0 and at most 1')
    }
    if (min_mac != null && (!(min_mac instanceof Number) || min_mac < 0)) {
        fail.call(analysis_id, 'regenie.min_mac', 'expected a non-negative number or null')
    }
    ['firth', 'firth_approx', 'firth_p_threshold'].each { option ->
        if (options.containsKey(option) && !methods.is_binary) {
            fail.call(analysis_id, "regenie.${option}", 'option is consumed by binary-trait REGENIE analyses only')
        }
    }
    if (options.containsKey('firth_p_threshold') && !firth) {
        fail.call(analysis_id, 'regenie.firth_p_threshold', "requires 'regenie.firth' to be true")
    }
    return [
        step1_bsize: step1_bsize.toInteger(),
        firth: firth,
        firth_approx: firth_approx,
        firth_p_threshold: firth_p_threshold,
        min_mac: min_mac,
    ]
}

def resolveGctaMethodOptions(analysis_id, options, methods, defaults, fail) {
    def capabilities = getMethodCapabilities()
    // The LD- and MAF-stratified plan settings live in this family because GCTA's own LD-score pass builds the
    // plan, but the plan itself is shared: one component plan now feeds the GCTA and the MPH stratified matrix
    // families. A row selecting only `mph_reml_ldms` therefore configures the plan through `gcta.*` and must
    // reach the per-option gates below rather than being turned away here for selecting no GCTA method.
    def stratified_heritability = getMethodTokensWithCapabilities([domain: 'heritability', component_model: 'ld_maf_stratified'])
    def selects_gcta = (methods.association_methods + methods.heritability_methods).any { method -> capabilities[method] && capabilities[method].option_family == 'gcta' }
    def selects_stratified = methods.heritability_methods.any { method -> method in stratified_heritability }
    if (options && !selects_gcta && !selects_stratified) {
        fail.call(analysis_id, "gcta.${options.keySet().first()}", 'analysis does not select a GCTA method')
    }
    def gcta_greml_estimators = getMethodTokensWithCapabilities([domain: 'heritability', option_family: 'gcta', estimator_family: 'reml'])
    def gcta_dense_heritability = getMethodTokensWithCapabilities([domain: 'heritability', option_family: 'gcta', input_backend: 'dense_grm'])
    def gcta_ldms_heritability = stratified_heritability
    def gcta_sparse_association = getMethodTokensWithCapabilities([domain: 'association', option_family: 'gcta', input_backend: 'sparse_grm'])
    def mph_matrix_heritability = getMethodTokensWithCapabilities([domain: 'heritability', option_family: 'mph'])
    if (options.containsKey('reml_no_constrain') && !methods.heritability_methods.any { method -> method in gcta_greml_estimators }) {
        fail.call(analysis_id, 'gcta.reml_no_constrain', "option is consumed by GCTA GREML estimators only, but this analysis selects neither 'gcta_greml' nor 'gcta_greml_ldms'")
    }
    ['grm_maf', 'grm_extract'].each { option ->
        if (options.containsKey(option) && !methods.heritability_methods.any { method -> method in gcta_dense_heritability }) {
            fail.call(analysis_id, "gcta.${option}", "option is consumed by 'gcta_greml' only, which this analysis does not select")
        }
    }
    // Both options restrict the GCTA dense predictor set and nothing else: matrix construction applies them to
    // `gcta_dense`/`gcta_sparse` only, and MPH's matrices are built over their own declared variant universe.
    // A row selecting both families would therefore publish a filtered GCTA estimate and an unfiltered MPH one
    // side by side under one analysis identifier with nothing in either native output saying so, which is the
    // failure mode a matched comparison exists to avoid. Refused rather than routed, following
    // `ldak.relatedness_filter`: making a `gcta.*` option govern MPH bytes would require renaming it out of
    // this family first.
    def selected_mph = methods.heritability_methods.findAll { method -> method in mph_matrix_heritability }
    ['grm_maf', 'grm_extract'].each { option ->
        if (options.containsKey(option) && selected_mph) {
            fail.call(analysis_id, "gcta.${option}", "option restricts the GCTA dense predictor set only; the MPH matrices this analysis also selects (${selected_mph.join(', ')}) are built over the declared MPH variant universe, so the two estimates would be published side by side over different predictor sets. Remove the option or drop the MPH selector(s) from this row")
        }
    }
    if (options.containsKey('sparse_cutoff') && !methods.association_methods.any { method -> method in gcta_sparse_association }) {
        fail.call(analysis_id, 'gcta.sparse_cutoff', "option is consumed by 'gcta_fastgwa' only, which this analysis does not select")
    }
    ['ld_score_region_kb', 'ld_bins', 'ldms_maf_edges'].each { option ->
        if (options.containsKey(option) && !methods.heritability_methods.any { method -> method in gcta_ldms_heritability }) {
            fail.call(analysis_id, "gcta.${option}", "option is consumed by the LD- and MAF-stratified estimators ${gcta_ldms_heritability.join(', ')} only, none of which this analysis selects")
        }
    }

    def grm_maf = options.containsKey('grm_maf') ? options.grm_maf : defaults.grm_maf
    def sparse_cutoff = options.containsKey('sparse_cutoff') ? options.sparse_cutoff : defaults.sparse_cutoff
    def ld_score_region_kb = options.containsKey('ld_score_region_kb') ? options.ld_score_region_kb : defaults.ld_score_region_kb
    def ld_bins = options.containsKey('ld_bins') ? options.ld_bins : defaults.ld_bins
    def ldms_maf_edges = options.containsKey('ldms_maf_edges') ? options.ldms_maf_edges : defaults.ldms_maf_edges
    def reml_no_constrain = options.containsKey('reml_no_constrain') ? options.reml_no_constrain : defaults.reml_no_constrain
    if (grm_maf != null && (!(grm_maf instanceof Number) || grm_maf < 0 || grm_maf > 0.5)) {
        fail.call(analysis_id, 'gcta.grm_maf', 'expected a number between 0 and 0.5 inclusive')
    }
    if (!(sparse_cutoff instanceof Number) || sparse_cutoff < 0 || sparse_cutoff > 1) {
        fail.call(analysis_id, 'gcta.sparse_cutoff', 'expected a number between 0 and 1 inclusive')
    }
    if (!(ld_score_region_kb instanceof Number) || ld_score_region_kb < 1 || ld_score_region_kb != ld_score_region_kb.toInteger()) {
        fail.call(analysis_id, 'gcta.ld_score_region_kb', 'expected one positive integer')
    }
    if (!(ld_bins instanceof Number) || ld_bins < 1 || ld_bins != ld_bins.toInteger()) {
        fail.call(analysis_id, 'gcta.ld_bins', 'expected one positive integer')
    }
    if (!(ldms_maf_edges instanceof List) || ldms_maf_edges.size() < 2 || !ldms_maf_edges.every { edge -> edge instanceof Number && edge >= 0 && edge <= 0.5 }) {
        fail.call(analysis_id, 'gcta.ldms_maf_edges', 'expected a numeric boundary list spanning 0 to 0.5')
    }
    if (ldms_maf_edges.first() != 0 || ldms_maf_edges.last() != 0.5) {
        fail.call(analysis_id, 'gcta.ldms_maf_edges', 'boundaries must start at 0 and end at 0.5')
    }
    if ((1..<ldms_maf_edges.size()).any { index -> ldms_maf_edges[index] <= ldms_maf_edges[index - 1] }) {
        fail.call(analysis_id, 'gcta.ldms_maf_edges', 'boundaries must be strictly increasing')
    }
    if (!(reml_no_constrain instanceof Boolean)) {
        fail.call(analysis_id, 'gcta.reml_no_constrain', 'expected a boolean')
    }
    return [
        grm_maf: grm_maf,
        grm_extract: resolveMethodResource(analysis_id, 'gcta', 'grm_extract', options, fail),
        reml_no_constrain: reml_no_constrain,
        sparse_cutoff: sparse_cutoff,
        ld_score_region_kb: ld_score_region_kb,
        ld_bins: ld_bins,
        ldms_maf_edges: ldms_maf_edges,
    ]
}

def resolveLdakMethodOptions(analysis_id, options, methods, defaults, fail) {
    def model = options.containsKey('model') ? options.model : defaults.model
    def power = options.containsKey('power') ? options.power : defaults.power
    def weights_policy = options.containsKey('weights_policy') ? options.weights_policy : defaults.weights_policy
    def relatedness_filter = options.containsKey('relatedness_filter') ? options.relatedness_filter : defaults.relatedness_filter
    def kvik_step1_subset = options.containsKey('kvik_step1_subset') ? options.kvik_step1_subset : defaults.kvik_step1_subset
    def kvik_step2_keep = resolveMethodResource(analysis_id, 'ldak', 'kvik_step2_keep', options, fail)
    if (!(model in ['human_default', 'custom'])) {
        fail.call(analysis_id, 'ldak.model', "expected 'human_default' or 'custom'")
    }
    if (!(power instanceof Number) || power < -2 || power > 0) {
        fail.call(analysis_id, 'ldak.power', 'expected a number between -2 and 0 inclusive')
    }
    if (model == 'human_default' && power != -0.25) {
        fail.call(analysis_id, 'ldak.power', "model 'human_default' fixes power at -0.25; set model to 'custom' to supply another power")
    }
    if (!(weights_policy in ['equal', 'provided'])) {
        fail.call(analysis_id, 'ldak.weights_policy', "expected 'equal' or 'provided'")
    }
    if (options.containsKey('weights') && weights_policy != 'provided') {
        fail.call(analysis_id, 'ldak.weights', "resource is only accepted when weights_policy is 'provided', got '${weights_policy}'")
    }
    if (weights_policy == 'provided' && !options.containsKey('weights')) {
        fail.call(analysis_id, 'ldak.weights', "weights_policy is 'provided' but no resource is supplied")
    }
    def weights = resolveMethodResource(analysis_id, 'ldak', 'weights', options, fail)
    if (!(relatedness_filter instanceof Boolean)) {
        fail.call(analysis_id, 'ldak.relatedness_filter', 'expected a boolean')
    }
    if (!(kvik_step1_subset in ['all', 'thin_common', 'provided'])) {
        fail.call(analysis_id, 'ldak.kvik_step1_subset', "expected 'all', 'thin_common' or 'provided'")
    }
    if (options.containsKey('predictor_extract') && kvik_step1_subset != 'provided') {
        fail.call(analysis_id, 'ldak.predictor_extract', "resource is only accepted when kvik_step1_subset is 'provided', got '${kvik_step1_subset}'")
    }
    if (kvik_step1_subset == 'provided' && !options.containsKey('predictor_extract')) {
        fail.call(analysis_id, 'ldak.predictor_extract', "kvik_step1_subset is 'provided' but no resource is supplied")
    }
    def predictor_extract = resolveMethodResource(analysis_id, 'ldak', 'predictor_extract', options, fail)

    def fast_repetitions = options.containsKey('fast_repetitions') ? options.fast_repetitions : defaults.fast_repetitions
    def fast_num_blocks = options.containsKey('fast_num_blocks') ? options.fast_num_blocks : defaults.fast_num_blocks
    def fast_seed = options.containsKey('fast_seed') ? options.fast_seed : defaults.fast_seed
    if (fast_repetitions != null && (!(fast_repetitions instanceof Number) || fast_repetitions < 1 || fast_repetitions != fast_repetitions.toInteger())) {
        fail.call(analysis_id, 'ldak.fast_repetitions', 'expected a positive integer or null')
    }
    if (fast_num_blocks != null && (!(fast_num_blocks instanceof Number) || fast_num_blocks != fast_num_blocks.toInteger())) {
        fail.call(analysis_id, 'ldak.fast_num_blocks', 'expected an integer or null')
    }
    if (fast_seed != null && (!(fast_seed instanceof Number) || fast_seed != fast_seed.toInteger())) {
        fail.call(analysis_id, 'ldak.fast_seed', 'expected an integer or null')
    }

    def ldak_kinship_heritability = getMethodTokensWithCapabilities([domain: 'heritability', input_backend: 'ldak_kinship'])
    def ldak_direct_heritability = getMethodTokensWithCapabilities([domain: 'heritability', option_family: 'ldak', input_backend: 'direct_plink1_genotypes'])
    def ldak_association = getMethodTokensWithCapabilities([domain: 'association', option_family: 'ldak'])
    def selects_ldak_kinship = methods.heritability_methods.any { method -> method in ldak_kinship_heritability }
    def selected_ldak_direct = methods.heritability_methods.findAll { method -> method in ldak_direct_heritability }
    def selects_ldak_kvik = methods.association_methods.any { method -> method in ldak_association }
    if (options && !selects_ldak_kinship && !selected_ldak_direct && !selects_ldak_kvik) {
        fail.call(analysis_id, "ldak.${options.keySet().first()}", 'analysis does not select an LDAK method')
    }
    // The model options describe how LDAK scales predictors, which both the kinship builder and the
    // direct-genotype estimators read; the unrelated-sample filter describes a kinship artifact only.
    if (options.containsKey('weights') && !selects_ldak_kinship && !selected_ldak_direct) {
        fail.call(analysis_id, 'ldak.weights', 'resource is consumed by LDAK kinship and direct-genotype heritability methods only, which this analysis does not select')
    }
    ['model', 'power', 'weights_policy'].each { option ->
        if (options.containsKey(option) && !selects_ldak_kinship && !selected_ldak_direct) {
            fail.call(analysis_id, "ldak.${option}", 'option is consumed by LDAK kinship and direct-genotype heritability methods only, which this analysis does not select')
        }
    }
    if (options.containsKey('relatedness_filter') && !selects_ldak_kinship) {
        fail.call(analysis_id, 'ldak.relatedness_filter', 'option is consumed by LDAK kinship methods only, which this analysis does not select')
    }
    ['fast_repetitions', 'fast_num_blocks', 'fast_seed'].each { option ->
        if (options.containsKey(option) && !selected_ldak_direct) {
            fail.call(analysis_id, "ldak.${option}", "option is consumed by the direct-genotype LDAK estimators ${ldak_direct_heritability.join(' and ')} only, which this analysis does not select")
        }
    }
    ['kvik_step1_subset', 'predictor_extract', 'kvik_step2_keep'].each { option ->
        if (options.containsKey(option) && !selects_ldak_kvik) {
            fail.call(analysis_id, "ldak.${option}", "option is consumed by 'ldak_kvik' only, which this analysis does not select")
        }
    }
    // A mixed row is allowed to keep the filter, because the kinship estimators genuinely consume it, but the
    // two families then estimate on different sample sets and nothing in either native output says so.
    if (relatedness_filter && selected_ldak_direct) {
        log.warn("[nf-core/gwas]: analysis '${analysis_id}' sets 'ldak.relatedness_filter' beside direct-genotype estimator(s) ${selected_ldak_direct.join(', ')}, which build no kinship and cannot consume an unrelated-sample keep list; the filter applies to the LDAK kinship estimators only, so the two estimates are computed on different sample sets")
    }
    return [
        model: model,
        power: power,
        weights_policy: weights_policy,
        weights: weights,
        relatedness_filter: relatedness_filter,
        kvik_step1_subset: kvik_step1_subset,
        predictor_extract: predictor_extract,
        kvik_step2_keep: kvik_step2_keep,
        fast_repetitions: fast_repetitions == null ? null : fast_repetitions.toInteger(),
        fast_num_blocks: fast_num_blocks == null ? null : fast_num_blocks.toInteger(),
        fast_seed: fast_seed == null ? null : fast_seed.toInteger(),
    ]
}

// MPH's controls split into two kinds and the resolver keeps them apart deliberately, because only one kind
// moves the answer. `seed`, `random_vectors` and `save_memory` govern the stochastic trace estimator: measured
// on the pinned image, the point estimate — not merely its standard error — moved across seeds, so
// a seed here changes the number that gets published. `iterations` and `tolerance` govern the deterministic solver.
//
// An unset `mph.seed` uses MPH's fixed integer 0, so an unseeded run is reproducible and its seed is recoverable
// from the log's option echo. `tolerance` is validated as any positive number rather than against MPH's own 1e-4 floor, because the
// floor is applied natively and this pipeline records what the researcher asked for beside what
// the tool did rather than pre-empting it.
def resolveMphMethodOptions(analysis_id, options, methods, defaults, fail) {
    def mph_heritability = getMethodTokensWithCapabilities([domain: 'heritability', option_family: 'mph'])
    if (options && !methods.heritability_methods.any { method -> method in mph_heritability }) {
        fail.call(analysis_id, "mph.${options.keySet().first()}", 'analysis does not select an MPH method')
    }
    return validateMphOptionValues(options, defaults) { option, reason -> fail.call(analysis_id, "mph.${option}", reason) }
}

// The MPH option value rules with no namespace and no selection knowledge, so the analysis namespace and the
// pair-request namespace apply one set of rules rather than two that can drift. The caller supplies a
// `fail(option_name, reason)` closure that renders its own namespace's diagnostic; the selection guard above
// stays in `resolveMphMethodOptions` because only the analysis namespace can declare an option for a method
// the row does not select. A pair request is keyed `<method>--<relationship_id>`, so its method is its identity.
def validateMphOptionValues(options, defaults, fail) {
    def iterations = options.containsKey('iterations') ? options.iterations : defaults.iterations
    def tolerance = options.containsKey('tolerance') ? options.tolerance : defaults.tolerance
    def random_vectors = options.containsKey('random_vectors') ? options.random_vectors : defaults.random_vectors
    def seed = options.containsKey('seed') ? options.seed : defaults.seed
    def save_memory = options.containsKey('save_memory') ? options.save_memory : defaults.save_memory
    if (iterations != null && (!(iterations instanceof Number) || iterations < 1 || iterations != iterations.toInteger())) {
        fail.call('iterations', 'expected a positive integer or null')
    }
    if (tolerance != null && (!(tolerance instanceof Number) || tolerance <= 0)) {
        fail.call('tolerance', 'expected a number greater than 0 or null')
    }
    if (random_vectors != null && (!(random_vectors instanceof Number) || random_vectors < 1 || random_vectors != random_vectors.toInteger())) {
        fail.call('random_vectors', 'expected a positive integer or null')
    }
    if (seed != null && (!(seed instanceof Number) || seed != seed.toInteger())) {
        fail.call('seed', 'expected an integer or null')
    }
    if (!(save_memory instanceof Boolean)) {
        fail.call('save_memory', 'expected a boolean')
    }
    return [
        iterations: iterations == null ? null : iterations.toInteger(),
        tolerance: tolerance,
        random_vectors: random_vectors == null ? null : random_vectors.toInteger(),
        seed: seed == null ? null : seed.toInteger(),
        save_memory: save_memory,
    ]
}

def readMethodOptionsDocument(method_options) {
    if (!method_options) {
        return [:]
    }
    def document_path = method_options.toString()
    def fail = { reason ->
        error("[nf-core/gwas] ERROR: Method-options document '${document_path}', analysis_id '<document>', option '<root>': ${reason}")
    }
    def document_file = file(method_options)
    if (!document_file.exists()) {
        fail.call('file does not exist')
    }

    def document = null
    try {
        document = new groovy.json.JsonSlurper().parseText(document_file.text)
    }
    catch (exception: Exception) {
        fail.call("malformed JSON (${exception.message})")
    }
    if (!(document instanceof Map)) {
        fail.call('expected an object keyed by analysis_id or by a supported request namespace')
    }
    return document
}

def getAnalysisOptionsDocument(method_options, document) {
    def namespaces = ['analyses', 'unary_requests', 'pair_requests', 'meta_requests']
    def uses_namespaces = document.keySet().any { key -> key in namespaces }
    if (!uses_namespaces) {
        return document
    }
    def unknown = document.keySet().findAll { key -> !(key in namespaces) }
    if (unknown) {
        error("[nf-core/gwas] ERROR: Method-options document '${method_options}' has unknown top-level namespace '${unknown.first()}'; accepted namespaces are ${namespaces.join(', ')}")
    }
    namespaces.each { namespace ->
        if (document.containsKey(namespace) && !(document[namespace] instanceof Map)) {
            error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}': expected an object")
        }
    }
    return document.analyses ?: [:]
}

// Parse and validate the analysis-owned part of the advanced method-options document.
def validateMethodOptions(method_options, analysis_rows, document = null) {
    def defaults = getMethodOptionDefaults()
    if (!method_options) {
        return analysis_rows.collectEntries { row -> [(row[0].id): defaults] }
    }

    document = document == null ? readMethodOptionsDocument(method_options) : document
    def document_path = method_options.toString()
    def fail = { analysis_id, option, reason ->
        error("[nf-core/gwas] ERROR: Method-options document '${document_path}', analysis_id '${analysis_id}', option '${option}': ${reason}")
    }
    def analysis_document = getAnalysisOptionsDocument(method_options, document)

    def analyses = analysis_rows.collectEntries { row ->
        def meta = row[0]
        [(meta.id): [
            association_methods: tokenizeMethodSelector(meta.association_methods),
            heritability_methods: tokenizeMethodSelector(meta.heritability_methods),
            is_binary: meta.trait_type == 'binary',
        ]]
    }
    def resolved = analyses.collectEntries { analysis_id, _methods -> [(analysis_id): defaults] }

    analysis_document.each { analysis_id, families ->
        if (!analyses.containsKey(analysis_id)) {
            fail.call(analysis_id, '<analysis>', 'analysis identifier is not declared in the analysis manifest')
        }
        if (!(families instanceof Map)) {
            fail.call(analysis_id, '<analysis>', 'expected a method-family object')
        }
        // Family acceptance, the option object check and per-family key validation all read the defaults map,
        // so a family's accepted names and its resolved shape come from one place. The resolver calls below
        // stay explicit: they are invoked in a fixed order because the first failure is the one a researcher
        // sees, and that attribution order is pinned by test.
        def unknown_families = families.keySet().findAll { family -> !defaults.containsKey(family) }
        if (unknown_families) {
            fail.call(analysis_id, unknown_families.first().toString(), "unknown method family; accepted families are ${defaults.keySet().join(', ')}")
        }

        // Object errors are attributed in declaration order and key errors in reverse, which is what makes the
        // first diagnostic a researcher sees stable across families rather than dependent on map iteration.
        def supplied = defaults.keySet().collectEntries { family -> [(family): families.containsKey(family) ? families[family] : [:]] }
        supplied.each { family, options ->
            if (!(options instanceof Map)) {
                fail.call(analysis_id, family, 'expected an option object')
            }
        }

        def operational = getMethodOptionOperationalKeys()
        defaults
            .keySet()
            .toList()
            .reverse()
            .each { family ->
                validateMethodOptionKeys(analysis_id, family, supplied[family], defaults[family].keySet().toList(), operational[family], fail)
            }

        def methods = analyses[analysis_id]
        def resolved_regenie = resolveRegenieMethodOptions(analysis_id, supplied.regenie, methods, defaults.regenie, fail)
        def resolved_gcta = resolveGctaMethodOptions(analysis_id, supplied.gcta, methods, defaults.gcta, fail)
        def resolved_ldak = resolveLdakMethodOptions(analysis_id, supplied.ldak, methods, defaults.ldak, fail)
        def resolved_mph = resolveMphMethodOptions(analysis_id, supplied.mph, methods, defaults.mph, fail)
        resolved[analysis_id] = [
            gcta: resolved_gcta,
            ldak: resolved_ldak,
            mph: resolved_mph,
            regenie: resolved_regenie,
        ]
    }

    return resolved
}
