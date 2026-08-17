// Validate linked cohort and analysis manifests as one pipeline input contract and emit canonical analyses.
// This workflow performs no tasks and therefore emits no versions.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// PLUGIN
include { samplesheetToList } from 'plugin/nf-schema'

workflow VALIDATE_GWAS_INPUT {
    take:
    cohort_manifest // channel: val(cohort_manifest)
    analysis_manifest // channel: val(analysis_manifest)
    method_options // channel: val(method_options)

    main:
    if (!cohort_manifest) {
        error("[nf-core/gwas] ERROR: --analysis_manifest '${analysis_manifest ?: ''}' was supplied but --cohort_manifest is missing; linked manifests require both")
    }
    if (!analysis_manifest) {
        error("[nf-core/gwas] ERROR: --cohort_manifest '${cohort_manifest}' was supplied but --analysis_manifest is missing; linked manifests require both")
    }

    def cohort_schema = "${projectDir}/assets/schema_cohort_manifest.json"
    def analysis_schema = "${projectDir}/assets/schema_analysis_manifest.json"
    validateSamplesheetHeader(cohort_manifest, cohort_schema, 'Cohort manifest')
    validateSamplesheetHeader(analysis_manifest, analysis_schema, 'Analysis manifest')

    def ch_analyses = channel.fromList(
        validateRelationalInput(
            samplesheetToList(cohort_manifest, cohort_schema),
            samplesheetToList(analysis_manifest, analysis_schema),
            cohort_manifest,
            analysis_manifest,
            cohort_schema,
            analysis_schema,
            method_options,
        )
    )

    emit:
    analyses = ch_analyses // channel: [ val(meta), [ path(genotype_file), ... ], path(phenotype), path(quant_covariates), path(cat_covariates), path(kvik_extract), path(ldak_weights) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// One registry owns selector domain, option family, matrix, prevalence and citation capabilities. Association
// entries also carry their GWASLab constructor mapping so the selector vocabulary cannot drift from it.
def getMethodRegistry() {
    return [
        plink2: [
            domain: 'association',
            citation_key: 'plink2',
            mapping: [common: [
                snpid: 'ID',
                chrom: '#CHROM',
                pos: 'POS',
                ea: 'A1',
                nea: 'REF',
                eaf: 'A1_FREQ',
                n: 'OBS_CT',
                beta: 'BETA',
                se: 'SE',
                p: 'P',
            ]],
        ],
        regenie: [
            domain: 'association',
            option_family: 'regenie',
            citation_key: 'regenie',
            mapping: [common: [
                snpid: 'ID',
                chrom: 'CHROM',
                pos: 'GENPOS',
                ea: 'ALLELE1',
                nea: 'ALLELE0',
                eaf: 'A1FREQ',
                n: 'N',
                beta: 'BETA',
                se: 'SE',
                mlog10p: 'LOG10P',
                readargs: [sep: ' '],
            ]],
        ],
        gcta_fastgwa: [
            domain: 'association',
            option_family: 'gcta',
            matrix_kind: 'gcta_sparse',
            citation_key: 'gcta_fastgwa',
            mapping: [common: [
                snpid: 'SNP',
                chrom: 'CHR',
                pos: 'POS',
                ea: 'A1',
                nea: 'A2',
                eaf: 'AF1',
                n: 'N',
                beta: 'BETA',
                se: 'SE',
                p: 'P',
            ]],
        ],
        ldak_kvik: [
            domain: 'association',
            option_family: 'ldak',
            citation_key: 'ldak_kvik',
            mapping: [
                common: [
                    snpid: 'Predictor',
                    chrom: 'Chromosome',
                    pos: 'Basepair',
                    ea: 'A1',
                    nea: 'A2',
                    eaf: 'EAF',
                    neff: 'N',
                    p: 'Wald_P',
                ],
                quantitative: [beta: 'Effect', se: 'SE'],
                binary: [beta: 'Approx_Log_OR', se: 'Approx_SE'],
            ],
        ],
        gcta_greml: [
            domain: 'heritability',
            option_family: 'gcta',
            matrix_kind: 'gcta_dense',
            consumes_population_prevalence: true,
            citation_key: 'gcta_greml',
        ],
        gcta_greml_ldms: [
            domain: 'heritability',
            option_family: 'gcta',
            matrix_kind: 'gcta_ldms',
            consumes_population_prevalence: true,
            citation_key: 'gcta_greml_ldms',
        ],
        ldak_reml: [
            domain: 'heritability',
            option_family: 'ldak',
            matrix_kind: 'ldak_kinship',
            consumes_population_prevalence: true,
            citation_key: 'ldak',
        ],
        ldak_he: [
            domain: 'heritability',
            option_family: 'ldak',
            matrix_kind: 'ldak_kinship',
            citation_key: 'ldak',
        ],
        ldak_pcgc: [
            domain: 'heritability',
            option_family: 'ldak',
            matrix_kind: 'ldak_kinship',
            consumes_population_prevalence: true,
            requires_population_prevalence: true,
            citation_key: 'ldak',
        ],
    ]
}

// Preserve the established public mapping shape while projecting it from the unified registry.
def getAssociationColumnMappings() {
    return getMethodRegistry()
        .findAll { _token, entry -> entry.domain == 'association' }
        .collectEntries { token, entry -> [(token): entry.mapping] }
}

def getAssociationColumnMappingJson(method, is_binary) {
    def entry = getAssociationColumnMappings()[method]
    if (!entry) {
        error("[nf-core/gwas] ERROR: no GWASLab column mapping is registered for association method '${method}'")
    }
    return groovy.json.JsonOutput.toJson(entry.common + (entry[is_binary ? 'binary' : 'quantitative'] ?: [:]))
}

// Project the stable capability interface without exposing the association mapping implementation.
def getMethodCapabilities() {
    return getMethodRegistry().collectEntries { token, entry ->
        [(token): entry.findAll { name, _value -> name != 'mapping' }]
    }
}

def getAssociationMethodTokens() {
    return getAssociationColumnMappings().keySet().toList()
}

def getHeritabilityMethodTokens() {
    return getMethodCapabilities()
        .findAll { _token, details -> details.domain == 'heritability' }
        .keySet()
        .toList()
}

def getGenotypeGroups() {
    return [
        plink2: ['pgen', 'psam', 'pvar'],
        plink1: ['bed', 'bim', 'fam'],
        vcf: ['vcf'],
    ]
}

// Derive nf-schema positional fields from each schema so adding a property cannot silently shift the tuple.
def getSamplesheetPositionalColumns(schema) {
    def properties = new groovy.json.JsonSlurper().parseText(file(schema).text).items.properties
    return properties.findAll { _column, definition -> !definition.containsKey('meta') }.keySet().toList()
}

// Optional columns remain mandatory headers. Diagnose missing, unexpected and repeated names before
// nf-schema can inject defaults or discard missing optional cells.
def validateSamplesheetHeader(samplesheet, schema, role) {
    def expected = new groovy.json.JsonSlurper().parseText(file(schema).text).items.properties.keySet().toList()
    def header_line = file(samplesheet).readLines().find { line -> line.trim() }
    def observed = header_line
        ? header_line.split(',', -1).collect { column -> column.trim().replaceAll(/^"|"$/, '') }
        : []
    def missing = expected.findAll { column -> !observed.contains(column) }
    def unexpected = observed.findAll { column -> !expected.contains(column) }.unique()
    def repeated = observed.countBy { column -> column }.findAll { _column, count -> count > 1 }.keySet().toList()

    def problems = []
    if (missing) {
        problems << "missing column headers: ${missing.collect { column -> "'${column}'" }.join(', ')}"
    }
    if (unexpected) {
        problems << "unexpected column headers: ${unexpected.collect { column -> "'${column}'" }.join(', ')}"
    }
    if (repeated) {
        problems << "repeated column headers: ${repeated.collect { column -> "'${column}'" }.join(', ')}"
    }
    if (problems) {
        error("[nf-core/gwas] ERROR: ${role} '${samplesheet}' header row 1 does not match the mandatory ${expected.size()}-column input contract.\n\n  - ${problems.join('\n  - ')}\n")
    }
}

// nf-schema represents an absent positional cell as an empty list; legal falsey values such as zero remain.
def normaliseCellValue(value) {
    if (value == null || (value instanceof Collection && value.isEmpty())) {
        return null
    }
    return value.toString().trim() ? value : null
}

def tokenizeMethodSelector(selector) {
    return selector ? selector.toString().tokenize(',').collect { token -> token.trim() }.findAll { token -> token } : []
}

def getMethodRoutes(association_methods, heritability_methods) {
    def capabilities = getMethodCapabilities()
    def selected = (association_methods + heritability_methods)
        .collect { method -> capabilities[method] }
        .findAll { details -> details }
    return [
        runs_heritability: selected.any { details -> details.domain == 'heritability' },
        consumes_population_prevalence: selected.any { details -> details.consumes_population_prevalence },
        runs_ldak_kvik: association_methods.contains('ldak_kvik'),
        runs_ldak_heritability: heritability_methods.any { method -> capabilities[method] && capabilities[method].option_family == 'ldak' },
        runs_ldak_pcgc: heritability_methods.any { method -> capabilities[method] && capabilities[method].requires_population_prevalence },
        runs_gcta: selected.any { details -> details.option_family == 'gcta' },
        runs_gcta_fastgwa: association_methods.contains('gcta_fastgwa'),
        runs_greml_ldms: heritability_methods.contains('gcta_greml_ldms'),
    ]
}

def getAnalysisSettings(meta) {
    return [
        population_prevalence: normaliseCellValue(meta.population_prevalence),
        case_value: normaliseCellValue(meta.case_value),
        control_value: normaliseCellValue(meta.control_value),
    ]
}

def validateMethodSelectors(association_methods, heritability_methods, reject) {
    [
        [column: 'association_methods', methods: association_methods, vocabulary: getAssociationMethodTokens()],
        [column: 'heritability_methods', methods: heritability_methods, vocabulary: getHeritabilityMethodTokens()],
    ].each { selector ->
        def unknown = selector.methods.findAll { method -> !selector.vocabulary.contains(method) }.unique()
        if (unknown) {
            reject.call(selector.column, "unknown method${unknown.size() > 1 ? 's' : ''} ${unknown.collect { method -> "'${method}'" }.join(', ')}, accepted values are ${selector.vocabulary.collect { method -> "'${method}'" }.join(', ')}")
        }
        def repeated = selector.methods.countBy { method -> method }.findAll { _method, count -> count > 1 }.keySet()
        if (repeated) {
            reject.call(selector.column, "method${repeated.size() > 1 ? 's' : ''} ${repeated.collect { method -> "'${method}'" }.join(', ')} listed more than once")
        }
    }
    if (!association_methods && !heritability_methods) {
        reject.call(['association_methods', 'heritability_methods'], 'row selects no method, populate one of them or remove the row')
    }
}

def validateGenotypeGroup(cells, reject) {
    def groups = getGenotypeGroups()
    def populated_groups = groups.findAll { _name, columns -> columns.any { column -> cells[column] } }
    if (!populated_groups) {
        reject.call(groups.values().flatten(), 'no genotype group is populated, supply exactly one of pgen/psam/pvar, bed/bim/fam or vcf')
    }
    else if (populated_groups.size() > 1) {
        populated_groups
            .keySet()
            .toList()
            .tail()
            .each { name ->
                def populated_column = groups[name].find { column -> cells[column] }
                reject.call(populated_column, 'a second genotype group is populated on this row, supply exactly one of pgen/psam/pvar, bed/bim/fam or vcf')
            }
    }
    else {
        def group_name = populated_groups.keySet().first()
        groups[group_name]
            .findAll { column -> !cells[column] }
            .each { column ->
                reject.call(column, "genotype group '${group_name}' is only partly populated, all of ${groups[group_name].join(', ')} are required together")
            }
    }
    return populated_groups.size() == 1 ? populated_groups.keySet().first() : null
}

def validateTraitColumns(is_binary, settings, reject) {
    if (is_binary) {
        if (settings.case_value == null) {
            reject.call('case_value', 'a binary trait must declare the value used for cases in the phenotype file')
        }
        if (settings.control_value == null) {
            reject.call('control_value', 'a binary trait must declare the value used for controls in the phenotype file')
        }
        if (settings.case_value != null && settings.control_value != null && settings.case_value.toString() == settings.control_value.toString()) {
            reject.call(['case_value', 'control_value'], 'binary case_value and control_value must be distinct source codes')
        }
    }
    else {
        if (settings.case_value != null) {
            reject.call('case_value', "case_value has no meaning on a quantitative trait, remove it or set trait_type to 'binary'")
        }
        if (settings.control_value != null) {
            reject.call('control_value', "control_value has no meaning on a quantitative trait, remove it or set trait_type to 'binary'")
        }
        if (settings.population_prevalence != null) {
            reject.call('population_prevalence', "prevalence has no meaning on a quantitative trait, remove it or set trait_type to 'binary'")
        }
    }
}

def validateMethodConditionedColumns(settings, routes, reject) {
    if (settings.population_prevalence != null && !routes.consumes_population_prevalence) {
        reject.call('population_prevalence', 'none of the selected estimators consumes it; select GCTA GREML, GCTA GREML-LDMS, LDAK REML or LDAK PCGC, or remove the prevalence')
    }
    if (settings.population_prevalence == null && routes.runs_ldak_pcgc) {
        reject.call('population_prevalence', "'ldak_pcgc' always estimates on the liability scale and requires a population prevalence")
    }
}

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
    if (options && !methods.association_methods.contains('regenie')) {
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
    def selected = (methods.association_methods + methods.heritability_methods).any { method -> capabilities[method] && capabilities[method].option_family == 'gcta' }
    if (options && !selected) {
        fail.call(analysis_id, "gcta.${options.keySet().first()}", 'analysis does not select a GCTA method')
    }
    if (options.containsKey('reml_no_constrain') && !methods.heritability_methods.any { method -> method in ['gcta_greml', 'gcta_greml_ldms'] }) {
        fail.call(analysis_id, 'gcta.reml_no_constrain', "option is consumed by GCTA GREML estimators only, but this analysis selects neither 'gcta_greml' nor 'gcta_greml_ldms'")
    }
    ['grm_maf', 'grm_extract'].each { option ->
        if (options.containsKey(option) && !methods.heritability_methods.contains('gcta_greml')) {
            fail.call(analysis_id, "gcta.${option}", "option is consumed by 'gcta_greml' only, which this analysis does not select")
        }
    }
    if (options.containsKey('sparse_cutoff') && !methods.association_methods.contains('gcta_fastgwa')) {
        fail.call(analysis_id, 'gcta.sparse_cutoff', "option is consumed by 'gcta_fastgwa' only, which this analysis does not select")
    }
    ['ld_score_region_kb', 'ld_bins', 'ldms_maf_edges'].each { option ->
        if (options.containsKey(option) && !methods.heritability_methods.contains('gcta_greml_ldms')) {
            fail.call(analysis_id, "gcta.${option}", "option is consumed by 'gcta_greml_ldms' only, which this analysis does not select")
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
    if (!(model in ['human_default', 'custom'])) {
        fail.call(analysis_id, 'ldak.model', "expected 'human_default' or 'custom'")
    }
    if (!(power instanceof Number) || power < -2 || power > 0) {
        fail.call(analysis_id, 'ldak.power', 'expected a number between -2 and 0 inclusive')
    }
    if (model == 'human_default' && power != -0.25) {
        fail.call(analysis_id, 'ldak.power', "model 'human_default' fixes power at -0.25; set model to 'custom' to supply another power")
    }
    if (!(weights_policy in ['equal', 'default', 'provided'])) {
        fail.call(analysis_id, 'ldak.weights_policy', "expected 'equal', 'default' or 'provided'")
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

    def capabilities = getMethodCapabilities()
    def selects_ldak_kinship = methods.heritability_methods.any { method -> capabilities[method] && capabilities[method].matrix_kind == 'ldak_kinship' }
    def selects_ldak_kvik = methods.association_methods.contains('ldak_kvik')
    if (options && !selects_ldak_kinship && !selects_ldak_kvik) {
        fail.call(analysis_id, "ldak.${options.keySet().first()}", 'analysis does not select an LDAK method')
    }
    if (options.containsKey('weights') && !selects_ldak_kinship) {
        fail.call(analysis_id, 'ldak.weights', 'resource is consumed by LDAK kinship methods only, which this analysis does not select')
    }
    ['model', 'power', 'weights_policy', 'relatedness_filter'].each { option ->
        if (options.containsKey(option) && !selects_ldak_kinship) {
            fail.call(analysis_id, "ldak.${option}", 'option is consumed by LDAK kinship methods only, which this analysis does not select')
        }
    }
    ['kvik_step1_subset', 'predictor_extract'].each { option ->
        if (options.containsKey(option) && !selects_ldak_kvik) {
            fail.call(analysis_id, "ldak.${option}", "option is consumed by 'ldak_kvik' only, which this analysis does not select")
        }
    }
    return [
        model: model,
        power: power,
        weights_policy: weights_policy,
        weights: weights,
        relatedness_filter: relatedness_filter,
        kvik_step1_subset: kvik_step1_subset,
        predictor_extract: predictor_extract,
    ]
}

// Parse and validate the advanced method-options document before any channels are constructed.
def validateMethodOptions(method_options, analysis_rows) {
    def defaults = getMethodOptionDefaults()
    if (!method_options) {
        return analysis_rows.collectEntries { row -> [(row[0].id): defaults] }
    }

    def document_path = method_options.toString()
    def fail = { analysis_id, option, reason ->
        error("[nf-core/gwas] ERROR: Method-options document '${document_path}', analysis_id '${analysis_id}', option '${option}': ${reason}")
    }
    def document_file = file(method_options)
    if (!document_file.exists()) {
        fail.call('<document>', '<root>', 'file does not exist')
    }

    def document = null
    try {
        document = new groovy.json.JsonSlurper().parseText(document_file.text)
    }
    catch (Exception exception) {
        fail.call('<document>', '<root>', "malformed JSON (${exception.message})")
    }
    if (!(document instanceof Map)) {
        fail.call('<document>', '<root>', 'expected an object keyed by analysis_id')
    }

    def analyses = analysis_rows.collectEntries { row ->
        def meta = row[0]
        [(meta.id): [
            association_methods: tokenizeMethodSelector(meta.association_methods),
            heritability_methods: tokenizeMethodSelector(meta.heritability_methods),
            is_binary: meta.trait_type == 'binary',
        ]]
    }
    def resolved = analyses.collectEntries { analysis_id, _methods -> [(analysis_id): defaults] }

    document.each { analysis_id, families ->
        if (!analyses.containsKey(analysis_id)) {
            fail.call(analysis_id, '<analysis>', 'analysis identifier is not declared in the analysis manifest')
        }
        if (!(families instanceof Map)) {
            fail.call(analysis_id, '<analysis>', 'expected a method-family object')
        }
        def unknown_families = families.keySet().findAll { family -> !defaults.containsKey(family) }
        if (unknown_families) {
            fail.call(analysis_id, unknown_families.first().toString(), 'unknown method family; accepted families are gcta, ldak and regenie')
        }

        def gcta = families.containsKey('gcta') ? families.gcta : [:]
        def ldak = families.containsKey('ldak') ? families.ldak : [:]
        def regenie = families.containsKey('regenie') ? families.regenie : [:]
        [[name: 'gcta', options: gcta], [name: 'ldak', options: ldak], [name: 'regenie', options: regenie]].each { family ->
            if (!(family.options instanceof Map)) {
                fail.call(analysis_id, family.name, 'expected an option object')
            }
        }

        validateMethodOptionKeys(analysis_id, 'regenie', regenie, defaults.regenie.keySet().toList(), ['step2_bsize', 'step1_mode', 'step1_jobs', 'lowmem'], fail)
        validateMethodOptionKeys(analysis_id, 'ldak', ldak, defaults.ldak.keySet().toList(), ['threads', 'jobs', 'partitions'], fail)
        validateMethodOptionKeys(analysis_id, 'gcta', gcta, defaults.gcta.keySet().toList(), ['gcta_grm_parts'], fail)

        def methods = analyses[analysis_id]
        def resolved_regenie = resolveRegenieMethodOptions(analysis_id, regenie, methods, defaults.regenie, fail)
        def resolved_gcta = resolveGctaMethodOptions(analysis_id, gcta, methods, defaults.gcta, fail)
        def resolved_ldak = resolveLdakMethodOptions(analysis_id, ldak, methods, defaults.ldak, fail)
        resolved[analysis_id] = [
            gcta: resolved_gcta,
            ldak: resolved_ldak,
            regenie: resolved_regenie,
        ]
    }
    return resolved
}

// Validate the lists as a linked contract and construct the canonical tuple consumed by GWAS.
def validateRelationalInput(cohort_rows, analysis_rows, cohort_manifest, analysis_manifest, cohort_schema, analysis_schema, method_options = null) {
    def cohort_columns = getSamplesheetPositionalColumns(cohort_schema)
    def analysis_columns = getSamplesheetPositionalColumns(analysis_schema)
    def errors = []
    def cohorts_by_id = [:]
    def line_by_analysis_id = [:]
    def validated_rows = []
    def method_options_by_analysis = validateMethodOptions(method_options, analysis_rows)

    cohort_rows.eachWithIndex { row, index ->
        def line = index + 2
        def cohort_meta = row[0]
        def cells = [cohort_columns, row[1..-1]].transpose().collectEntries()
        def cohort_id = cohort_meta.cohort
        def reject = { field, message ->
            def named = field instanceof List ? field : [field]
            def label = named.size() > 1
                ? "fields ${named.collect { name -> "'${name}'" }.join(', ')}"
                : "field '${named.first()}'"
            errors << "  - ${cohort_manifest} row ${line} (cohort_id '${cohort_id}'), ${label}: ${message}"
        }

        def genotype_format = validateGenotypeGroup(cells, reject)
        def genotype_files = genotype_format
            ? getGenotypeGroups()[genotype_format].collect { column -> cells[column] }
            : []
        def definition = [
            genome_build: cohort_meta.build,
            ancestry: cohort_meta.ancestry,
        ] + getGenotypeGroups().values().flatten().collectEntries { field -> [(field): normaliseCellValue(cells[field])?.toString() ?: ''] }
        def known = cohorts_by_id[cohort_id]

        if (known) {
            if (known.definition == definition) {
                reject.call('cohort_id', "duplicate cohort_id '${cohort_id}', identical definition on row ${known.line}")
            }
            else {
                def differing_fields = definition.keySet().findAll { field -> known.definition[field] != definition[field] }
                reject.call('cohort_id', "duplicate cohort_id '${cohort_id}' conflicts with row ${known.line} in field${differing_fields.size() > 1 ? 's' : ''} ${differing_fields.collect { field -> "'${field}'" }.join(', ')}")
            }
        }
        else {
            cohorts_by_id[cohort_id] = [
                line: line,
                meta: cohort_meta,
                genotype_format: genotype_format,
                genotype_files: genotype_files,
                definition: definition,
            ]
        }
    }

    analysis_rows.eachWithIndex { row, index ->
        def line = index + 2
        def analysis_meta = row[0]
        def cells = [analysis_columns, row[1..-1]].transpose().collectEntries()
        def analysis_id = analysis_meta.id
        def cohort_id = analysis_meta.cohort
        def reject = { field, message ->
            def named = field instanceof List ? field : [field]
            def label = named.size() > 1
                ? "fields ${named.collect { name -> "'${name}'" }.join(', ')}"
                : "field '${named.first()}'"
            errors << "  - ${analysis_manifest} row ${line} (analysis_id '${analysis_id}'), ${label}: ${message}"
        }

        if (line_by_analysis_id.containsKey(analysis_id)) {
            reject.call('analysis_id', "duplicate analysis_id '${analysis_id}', already declared on row ${line_by_analysis_id[analysis_id]}")
        }
        else {
            line_by_analysis_id[analysis_id] = line
        }

        def cohort = cohorts_by_id[cohort_id]
        if (!cohort) {
            reject.call('cohort_id', "undefined cohort_id '${cohort_id}', not declared in cohort manifest '${cohort_manifest}'")
        }

        def association_methods = tokenizeMethodSelector(analysis_meta.association_methods)
        def heritability_methods = tokenizeMethodSelector(analysis_meta.heritability_methods)
        validateMethodSelectors(association_methods, heritability_methods, reject)
        def routes = getMethodRoutes(association_methods, heritability_methods)
        def settings = getAnalysisSettings(analysis_meta)
        def is_binary = analysis_meta.trait_type == 'binary'
        validateTraitColumns(is_binary, settings, reject)
        validateMethodConditionedColumns(settings, routes, reject)

        if (cohort) {
            validated_rows << [
                analysis_meta + [
                    build: cohort.meta.build,
                    ancestry: cohort.meta.ancestry,
                    association_methods: association_methods,
                    heritability_methods: heritability_methods,
                    genotype_format: cohort.genotype_format,
                    is_binary: is_binary,
                    has_covariates: cells.quant_covariates || cells.cat_covariates ? true : false,
                    case_value: settings.case_value == null ? null : settings.case_value.toString(),
                    control_value: settings.control_value == null ? null : settings.control_value.toString(),
                    population_prevalence: settings.population_prevalence,
                    method_options: method_options_by_analysis[analysis_id],
                ],
                cohort.genotype_files,
                cells.phenotype,
                cells.quant_covariates,
                cells.cat_covariates,
                method_options_by_analysis[analysis_id].ldak.predictor_extract,
                method_options_by_analysis[analysis_id].ldak.weights,
            ]
        }
    }

    if (errors) {
        error("[nf-core/gwas] ERROR: Validation of linked manifests failed!\n\n${errors.join('\n')}\n")
    }
    return validated_rows
}
