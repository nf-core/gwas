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
        validateSamplesheetHeader(cohort_manifest, cohort_schema, 'Cohort manifest')
        validateSamplesheetHeader(analysis_manifest, analysis_schema, 'Analysis manifest')
    }
    if (summary_statistics_manifest) {
        validateSamplesheetHeader(summary_statistics_manifest, summary_statistics_schema, 'Summary-statistics manifest')
    }
    if (relationship_manifest) {
        validateSamplesheetHeader(relationship_manifest, relationship_schema, 'Relationship manifest')
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

    emit:
    analyses           = ch_analyses // channel: [ val(meta), [ path(genotype_file), ... ], path(phenotype), path(quant_covariates), path(cat_covariates), path(kvik_extract), path(ldak_weights) ]
    summary_statistics = ch_summary_statistics // channel: [ val(meta), path(source) ]
    relationships      = ch_relationships // channel: [ val(meta), [ path(genotype_file), ... ], path(pair_quant_covariates), path(pair_cat_covariates) ]
    unary_requests     = ch_unary_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ]
    pair_requests      = ch_pair_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ]
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
        gcta_bivariate_reml: [
            domain: 'pairwise',
            endpoint_domain: 'analysis',
            option_family: 'gcta',
            matrix_kind: 'gcta_dense',
            consumes_population_prevalence: true,
            citation_key: 'gcta_bivariate_reml',
        ],
        gcta_bivariate_reml_ldms: [
            domain: 'pairwise',
            endpoint_domain: 'analysis',
            option_family: 'gcta',
            matrix_kind: 'gcta_ldms',
            consumes_population_prevalence: true,
            citation_key: 'gcta_bivariate_reml',
            citation_keys: ['gcta_bivariate_reml', 'gcta_greml_ldms'],
        ],
        ldak_sumher: [
            domain: 'summary_unary',
            endpoint_domain: 'summary_statistics',
            option_family: 'ldak',
            reference_family: 'ldak',
            consumes_population_prevalence: true,
            consumes_sample_prevalence: true,
            citation_key: 'ldak_sumstats',
        ],
        ldak_sumcors: [
            domain: 'pairwise',
            endpoint_domain: 'summary_statistics',
            option_family: 'ldak',
            reference_family: 'ldak',
            consumes_population_prevalence: true,
            consumes_sample_prevalence: true,
            citation_key: 'ldak_sumstats',
        ],
        ldsc_h2: [
            domain: 'summary_unary',
            endpoint_domain: 'summary_statistics',
            option_family: 'ldsc',
            reference_family: 'ldsc',
            consumes_population_prevalence: true,
            consumes_sample_prevalence: true,
            citation_key: 'ldsc',
        ],
        ldsc_rg: [
            domain: 'pairwise',
            endpoint_domain: 'summary_statistics',
            option_family: 'ldsc',
            reference_family: 'ldsc',
            consumes_population_prevalence: true,
            consumes_sample_prevalence: true,
            citation_key: 'ldsc',
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

def getSummaryUnaryMethodTokens() {
    return getMethodCapabilities()
        .findAll { _token, details -> details.domain == 'summary_unary' }
        .keySet()
        .toList()
}

def getRelationshipMethodTokens() {
    return getMethodCapabilities()
        .findAll { _token, details -> details.domain == 'pairwise' }
        .keySet()
        .toList()
}

def getSummaryStatisticsId(analysis_id, association_method) {
    return "${analysis_id}--${association_method}".toString()
}

def getInternalSummaryMetadata(meta, association_method) {
    def summary_statistics_id = getSummaryStatisticsId(meta.id, association_method)
    return [
        id: summary_statistics_id,
        summary_statistics_id: summary_statistics_id,
        trait: meta.trait,
        trait_id: meta.trait,
        trait_type: meta.trait_type,
        is_binary: meta.is_binary,
        population_prevalence: meta.population_prevalence,
        sample_prevalence: meta.sample_prevalence,
        build: meta.build,
        ancestry: meta.ancestry,
        source_kind: 'pipeline_generated',
        source_mode: 'raw',
        source_format: "pipeline_${association_method}",
        source_method: association_method,
        source_release: null,
        source_name: "${meta.id}.${association_method}",
        producer_analysis_id: meta.id,
        producer_association_method: association_method,
        canonical_contract_version: 'nfcore_gwas_canonical_v1',
        transformation: 'gwaslab_harmonised',
        harmonization: [tool: 'GWASLab', version: '4.1.9'],
        access_constraints: null,
    ]
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
        sample_prevalence: normaliseCellValue(meta.sample_prevalence),
        case_value: normaliseCellValue(meta.case_value),
        control_value: normaliseCellValue(meta.control_value),
    ]
}

def validateMethodSelectors(association_methods, heritability_methods, reject, allow_empty = false) {
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
    if (!allow_empty && !association_methods && !heritability_methods) {
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
        if (settings.sample_prevalence != null) {
            reject.call('sample_prevalence', "sample prevalence has no meaning on a quantitative trait, remove it or set trait_type to 'binary'")
        }
    }
}

def validateMethodConditionedColumns(settings, routes, reject, downstream_consumes_population_prevalence = false) {
    if (settings.population_prevalence != null && !routes.consumes_population_prevalence && !downstream_consumes_population_prevalence) {
        reject.call('population_prevalence', 'none of the selected estimators consumes it; select a liability-aware individual, summary or pair method, or remove the prevalence')
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
    catch (Exception exception) {
        fail.call("malformed JSON (${exception.message})")
    }
    if (!(document instanceof Map)) {
        fail.call('expected an object keyed by analysis_id or by a supported request namespace')
    }
    return document
}

def getAnalysisOptionsDocument(method_options, document) {
    def namespaces = ['analyses', 'unary_requests', 'pair_requests']
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

def readReferenceCatalog(reference_catalog) {
    if (!reference_catalog) {
        return [:]
    }
    def catalog_path = reference_catalog.toString()
    def fail = { bundle_id, field, reason ->
        error("[nf-core/gwas] ERROR: Reference catalog '${catalog_path}', reference_bundle_id '${bundle_id}', field '${field}': ${reason}")
    }
    def catalog_file = file(reference_catalog)
    if (!catalog_file.exists()) {
        fail.call('<document>', '<root>', 'file does not exist')
    }
    def document = null
    try {
        document = new groovy.json.JsonSlurper().parseText(catalog_file.text)
    }
    catch (Exception exception) {
        fail.call('<document>', '<root>', "malformed JSON (${exception.message})")
    }
    if (!(document instanceof Map)) {
        fail.call('<document>', '<root>', 'expected an object with optional ldsc and ldak family objects')
    }
    def unknown_families = document.keySet().findAll { family -> !(family in ['ldsc', 'ldak']) }
    if (unknown_families) {
        fail.call('<document>', unknown_families.first().toString(), "unknown family; accepted families are 'ldsc' and 'ldak'")
    }

    def resolved = [:]
    ['ldsc', 'ldak'].each { family ->
        def bundles = document[family] ?: [:]
        if (!(bundles instanceof Map)) {
            fail.call('<document>', family, 'expected an object keyed by reference_bundle_id')
        }
        bundles.each { bundle_id, definition ->
            if (!(bundle_id instanceof String) || !(bundle_id ==~ /^\S+$/)) {
                fail.call(bundle_id, '<id>', 'expected one non-empty identifier without spaces')
            }
            if (resolved.containsKey(bundle_id)) {
                fail.call(bundle_id, '<id>', "identifier is already declared in family '${resolved[bundle_id].family}'")
            }
            if (!(definition instanceof Map)) {
                fail.call(bundle_id, '<bundle>', 'expected an object')
            }
            def required = family == 'ldsc'
                ? ['genome_build', 'ancestry', 'variant_id_system', 'hapmap3_snplist', 'reference_ld_scores', 'regression_weights']
                : ['genome_build', 'ancestry', 'variant_id_system', 'model', 'tagging_file']
            def optional = family == 'ldsc'
                ? ['hapmap3_sha256', 'reference_ld_scores_sha256', 'regression_weights_sha256']
                : ['tagging_sha256']
            def unknown = definition.keySet().findAll { field -> !(field in required + optional) }
            if (unknown) {
                fail.call(bundle_id, unknown.first().toString(), "unknown field; accepted fields are ${(required + optional).join(', ')}")
            }
            required.each { field ->
                if (!definition.containsKey(field) || definition[field] == null || !definition[field].toString().trim()) {
                    fail.call(bundle_id, field, 'required value is missing')
                }
            }
            if (!(definition.genome_build in ['GRCh37', 'GRCh38'])) {
                fail.call(bundle_id, 'genome_build', "expected 'GRCh37' or 'GRCh38'")
            }
            if (family == 'ldak' && !(definition.model in ['BLD-LDAK', 'Baseline-LD-v2.2', 'LDAK-Thin', 'Uniform-GCTA', 'Human-Default'])) {
                fail.call(bundle_id, 'model', "unsupported first-release model; expected 'BLD-LDAK', 'Baseline-LD-v2.2', 'LDAK-Thin', 'Uniform-GCTA' or 'Human-Default'")
            }
            optional
                .findAll { field -> definition.containsKey(field) }
                .each { field ->
                    if (!(definition[field].toString() ==~ /^[a-fA-F0-9]{64}$/)) {
                        fail.call(bundle_id, field, 'expected one SHA-256 digest containing exactly 64 hexadecimal characters')
                    }
                }

            def role_fields = family == 'ldsc'
                ? ['hapmap3_snplist', 'reference_ld_scores', 'regression_weights']
                : ['tagging_file']
            def resources = role_fields.collectEntries { field ->
                def resource = file(definition[field])
                if (!resource.exists()) {
                    fail.call(bundle_id, field, "resource path '${definition[field]}' does not exist")
                }
                if (field in ['hapmap3_snplist', 'tagging_file'] && !resource.toFile().isFile()) {
                    fail.call(bundle_id, field, "resource path '${definition[field]}' must be a file")
                }
                [(field): resource]
            }
            resolved[bundle_id] = [
                id: bundle_id,
                family: family,
                genome_build: definition.genome_build,
                ancestry: definition.ancestry,
                variant_id_system: definition.variant_id_system,
                model: family == 'ldak' ? definition.model : null,
                role_names: role_fields.collectEntries { field -> [(field): resources[field].name] },
                declared_checksums: optional.collectEntries { field -> [(field): definition[field] ?: null] },
                resources: resources,
            ]
        }
    }
    return resolved
}

def getLdscProtectedNativeOptionMatches(method, option_names) {
    def protected_by_method = [
        ldsc_h2: [
            wrapper_owned: ['--h2', '--ref-ld-chr', '--w-ld-chr', '--samp-prev', '--pop-prev', '--out'],
            typed_resource: [
                '--annot',
                '--bfile',
                '--cts-bin',
                '--extract',
                '--frqfile',
                '--frqfile-chr',
                '--h2-cts',
                '--keep',
                '--print-snps',
                '--ref-ld',
                '--ref-ld-chr-cts',
                '--w-ld',
            ],
            alternate_operation: ['--l2', '--rg'],
        ],
        ldsc_rg: [
            wrapper_owned: ['--rg', '--ref-ld-chr', '--w-ld-chr', '--samp-prev', '--pop-prev', '--out'],
            typed_resource: [
                '--annot',
                '--bfile',
                '--cts-bin',
                '--extract',
                '--frqfile',
                '--frqfile-chr',
                '--h2-cts',
                '--keep',
                '--print-snps',
                '--ref-ld',
                '--ref-ld-chr-cts',
                '--w-ld',
            ],
            alternate_operation: ['--l2', '--h2'],
        ],
    ]
    def protection_sets = protected_by_method[method] ?: [:]
    return option_names.collectEntries { option_name ->
        def exact = protection_sets.collectMany { kind, options -> options.findAll { protected_option -> protected_option == option_name }.collect { protected_option -> [kind: kind, option: protected_option] } }
        def prefix = option_name.startsWith('--') && option_name.size() > 2
            ? protection_sets.collectMany { kind, options ->
                options
                    .findAll { protected_option -> protected_option.startsWith(option_name) }
                    .collect { protected_option -> [kind: kind, option: protected_option] }
            }
            : []
        def matches = exact ?: prefix
        [(option_name): matches ? [exact: !exact.isEmpty(), matches: matches.unique()] : null]
    }
}

def validateSummaryNativeArgumentTokens(method_options, namespace, request_id, method, native_args) {
    def fail = { reason ->
        error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request_id}', option 'native_args': ${reason}")
    }
    if (!(native_args instanceof List)) {
        fail.call('expected an array of individual command-line argument tokens')
    }
    if (native_args && !native_args.first().toString().startsWith('--')) {
        fail.call("the first token must be a native option beginning with '--'")
    }
    def reserved_by_method = [
        ldak_sumher: [
            '--sum-hers',
            '--sum-cors',
            '--summary',
            '--summary2',
            '--tagfile',
            '--out',
            '--threads',
            '--max-threads',
            '--prevalence',
            '--ascertainment',
            '--prevalence2',
            '--ascertainment2',
        ],
        ldak_sumcors: [
            '--sum-hers',
            '--sum-cors',
            '--summary',
            '--summary2',
            '--tagfile',
            '--out',
            '--threads',
            '--max-threads',
            '--prevalence',
            '--ascertainment',
            '--prevalence2',
            '--ascertainment2',
        ],
    ]
    def undeclared_file_options = [
        ldak_sumher: [
            '--alternative-tags',
            '--categories',
            '--exclude',
            '--extract',
            '--keep',
            '--labels',
            '--matrix',
            '--remove',
            '--weights',
        ],
        ldak_sumcors: [
            '--alternative-tags',
            '--categories',
            '--exclude',
            '--extract',
            '--keep',
            '--labels',
            '--matrix',
            '--remove',
            '--weights',
        ],
    ]
    native_args.eachWithIndex { token, index ->
        if (!(token instanceof String) || !token) {
            fail.call("token ${index + 1} must be a non-empty string")
        }
        if (!(token ==~ /^[A-Za-z0-9_.:+,@%=-]+$/)) {
            fail.call("token ${index + 1} '${token}' contains whitespace, shell syntax or a path separator; pass individual non-file native tokens only")
        }
        if (token ==~ /^[A-Za-z_][A-Za-z0-9_]*=.*/) {
            fail.call("token ${index + 1} '${token}' resembles an environment assignment; native arguments cannot alter the task environment")
        }
        def option_name = token.contains('=') ? token.substring(0, token.indexOf('=')) : token
        if (method in ['ldsc_h2', 'ldsc_rg']) {
            def protection = getLdscProtectedNativeOptionMatches(method, [option_name])[option_name]
            if (protection) {
                if (!protection.exact) {
                    def matched_options = protection.matches.collect { match -> match.option }.unique().sort().join(', ')
                    fail.call("token ${index + 1} '${token}' is a protected LDSC option abbreviation matching ${matched_options}; abbreviated options cannot bypass wrapper-owned invocation, typed-resource or primary-operation controls")
                }
                def kind = protection.matches.first().kind
                if (kind == 'wrapper_owned') {
                    fail.call("token ${index + 1} '${token}' conflicts with wrapper-owned invocation mechanics")
                }
                if (kind == 'typed_resource') {
                    fail.call("token ${index + 1} '${token}' requires a typed staged resource, but this request architecture declares no such file role")
                }
                fail.call("token ${index + 1} '${token}' selects a different primary operation from wrapper-owned method '${method}'")
            }
        }
        else {
            if (option_name in (reserved_by_method[method] ?: [])) {
                fail.call("token ${index + 1} '${token}' conflicts with wrapper-owned invocation mechanics")
            }
            if (option_name in (undeclared_file_options[method] ?: [])) {
                fail.call("token ${index + 1} '${token}' requires a typed staged resource, but this request architecture declares no such file role")
            }
        }
        def argument_value = token.contains('=') ? token.substring(token.indexOf('=') + 1) : token
        if (!token.startsWith('--') || token.contains('=')) {
            def candidate = file(argument_value)
            def looks_like_file = argument_value ==~ /(?i).+\.(gz|bgz|txt|tsv|csv|list|tagging|annot|l2\.ldscore|weights|sumstats)/
            if (candidate.exists() || looks_like_file) {
                fail.call("token ${index + 1} '${token}' resembles an undeclared file input; file-taking native options require a typed staged resource")
            }
        }
    }
    if (method in ['ldak_sumher', 'ldak_sumcors']) {
        def option_names = native_args
            .findAll { token -> token instanceof String && token.startsWith('--') }
            .collect { token -> token.split('=', 2)[0] }
        if (option_names.contains('--cutoff') && option_names.contains('--truncate')) {
            fail.call("'--cutoff' and '--truncate' are mutually exclusive LDAK large-effect policies")
        }
    }
    return native_args
}

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
        if (request.meta.matrix_kind == 'gcta_ldms') {
            accepted += ['ld_score_region_kb', 'ld_bins', 'ldms_maf_edges']
        }
        def unknown = options.keySet().findAll { option -> !(option in accepted) }
        if (unknown) {
            error("[nf-core/gwas] ERROR: Method-options document '${method_options}', namespace '${namespace}', request_id '${request.request_id}', option '${unknown.first()}': unknown option; accepted options are ${accepted.join(', ')}")
        }
        def native_args = request.reference_family
            ? validateSummaryNativeArgumentTokens(method_options, namespace, request.request_id, request.method, options.native_args ?: [])
            : validateNativeArgumentTokens(method_options, request.request_id, options.native_args ?: [])
        def matrix_settings = request.meta.matrix_settings ?: [:]
        if (request.meta.matrix_kind == 'gcta_ldms') {
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
            if (request.method == 'ldak_sumcors' && !(bundle.model in ['LDAK-Thin', 'Uniform-GCTA', 'Human-Default'])) {
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
            matrix_settings: matrix_settings,
            reference_bundle_id: bundle ? bundle.id : null,
            reference_family: bundle ? bundle.family : null,
            reference_metadata: bundle
                ? bundle.subMap(['id', 'family', 'genome_build', 'ancestry', 'variant_id_system', 'model', 'role_names', 'declared_checksums'])
                : null,
            reference_validation: bundle ? 'structural_availability_only' : null,
        ]
        resolved << [request: request, meta: resolved_meta, bundle: bundle]
    }
    return resolved
}

def validateNativeArgumentTokens(method_options, request_id, native_args) {
    def fail = { reason ->
        error("[nf-core/gwas] ERROR: Method-options document '${method_options}', pair_request_id '${request_id}', option 'native_args': ${reason}")
    }
    if (!(native_args instanceof List)) {
        fail.call('expected an array of individual command-line argument tokens')
    }
    if (native_args && !native_args.first().toString().startsWith('--')) {
        fail.call("the first token must be a native option beginning with '--'")
    }

    def reserved = [
        '--reml-bivar',
        '--reml-bivar-prevalence',
        '--bfile',
        '--pfile',
        '--mbfile',
        '--mpfile',
        '--grm',
        '--mgrm',
        '--grm-gz',
        '--mgrm-gz',
        '--pheno',
        '--mpheno',
        '--qcovar',
        '--covar',
        '--keep',
        '--remove',
        '--extract',
        '--exclude',
        '--update-sex',
        '--dosage-mach',
        '--dosage-beagle',
        '--raw-files',
        '--out',
        '--thread-num',
        '--prevalence',
        '--help',
        '--version',
    ]
    def primary_operations = [
        '--reml',
        '--reml-ldms',
        '--make-grm',
        '--make-grm-part',
        '--make-bK-sparse',
        '--fastGWA-mlm',
        '--fastGWA-mlm-binary',
        '--mlma',
        '--cojo-slct',
        '--cojo-joint',
        '--HEreg',
        '--pca',
        '--pca-loading',
        '--freq',
        '--ld',
        '--ld-score',
        '--sblup',
        '--blup-snp',
        '--simu-qt',
        '--simu-cc',
        '--simu-causal-loci',
        '--GTDT',
    ]

    native_args.eachWithIndex { token, index ->
        if (!(token instanceof String) || !token) {
            fail.call("token ${index + 1} must be a non-empty string")
        }
        if (!(token ==~ /^[A-Za-z0-9_.:+,@%=-]+$/)) {
            fail.call("token ${index + 1} '${token}' contains whitespace, shell syntax or a path separator; pass individual non-file native tokens only")
        }
        if (token ==~ /^[A-Za-z_][A-Za-z0-9_]*=.*/) {
            fail.call("token ${index + 1} '${token}' resembles an environment assignment; native arguments cannot alter the task environment")
        }
        def option_name = token.contains('=') ? token.substring(0, token.indexOf('=')) : token
        if (option_name in reserved) {
            fail.call("token ${index + 1} '${token}' conflicts with wrapper-owned invocation mechanics")
        }
        if (option_name in primary_operations || option_name.startsWith('--make-') || option_name.startsWith('--fastGWA') || option_name.startsWith('--mlma') || option_name.startsWith('--cojo-') || option_name.startsWith('--simu-')) {
            fail.call("token ${index + 1} '${token}' selects a different primary GCTA operation")
        }
        def argument_value = token.contains('=') ? token.substring(token.indexOf('=') + 1) : token
        if (!token.startsWith('--') || token.contains('=')) {
            def candidate = file(argument_value)
            def looks_like_file = argument_value ==~ /(?i).+\.(bed|bim|fam|pgen|pvar|psam|vcf|bcf|gz|bgz|txt|tsv|csv|list|keep|remove|grm|mgrm|phen|pheno|covar|qcovar|dat)/
            if (candidate.exists() || looks_like_file) {
                fail.call("token ${index + 1} '${token}' resembles an undeclared file input; file-taking native options require a typed staged resource")
            }
        }
    }
    return native_args
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

// Validate the lists as a linked contract and construct the canonical tuple consumed by GWAS.
def validateRelationalInput(cohort_rows, analysis_rows, summary_statistics_rows, relationship_rows, cohort_manifest, analysis_manifest, summary_statistics_manifest, relationship_manifest, cohort_schema, analysis_schema, summary_statistics_schema, relationship_schema, reference_catalog, method_options = null) {
    def cohort_columns = getSamplesheetPositionalColumns(cohort_schema)
    def analysis_columns = getSamplesheetPositionalColumns(analysis_schema)
    def summary_statistics_columns = getSamplesheetPositionalColumns(summary_statistics_schema)
    def relationship_columns = getSamplesheetPositionalColumns(relationship_schema)
    def errors = []
    def cohorts_by_id = [:]
    def analyses_by_id = [:]
    def summaries_by_id = [:]
    def line_by_analysis_id = [:]
    def line_by_summary_statistics_id = [:]
    def line_by_relationship_id = [:]
    def validated_analyses = []
    def validated_external_summaries = []
    def declared_summaries = []
    def generated_summaries_by_id = [:]
    def gcta_primary_requests = []
    def summary_unary_primary_requests = []
    def summary_pair_primary_requests = []
    def options_document = readMethodOptionsDocument(method_options)
    def method_options_by_analysis = validateMethodOptions(method_options, analysis_rows, options_document)
    def reference_bundles = readReferenceCatalog(reference_catalog)
    def pair_prevalence_analysis_ids = relationship_rows
        .findAll { row ->
            tokenizeMethodSelector(row[0].relationship_methods).any { method ->
                def capability = getMethodCapabilities()[method]
                capability && capability.domain == 'pairwise' && capability.consumes_population_prevalence
            }
        }
        .collectMany { row -> [normaliseCellValue(row[0].left_analysis_id), normaliseCellValue(row[0].right_analysis_id)] }
        .findAll { analysis_id -> analysis_id }
        .collect { analysis_id -> analysis_id.toString() } as Set
    def pair_prevalence_summary_statistics_ids = relationship_rows
        .findAll { row ->
            tokenizeMethodSelector(row[0].relationship_methods).any { method ->
                def capability = getMethodCapabilities()[method]
                capability && capability.domain == 'pairwise' && capability.consumes_population_prevalence
            }
        }
        .collectMany { row -> [normaliseCellValue(row[0].left_summary_statistics_id), normaliseCellValue(row[0].right_summary_statistics_id)] }
        .findAll { summary_statistics_id -> summary_statistics_id }
        .collect { summary_statistics_id -> summary_statistics_id.toString() } as Set
    def summary_prevalence_analysis_ids = summary_statistics_rows
        .findAll { row ->
            def summary_statistics_id = normaliseCellValue(row[0].id)?.toString()
            normaliseCellValue(row[0].producer_analysis_id) && (
                tokenizeMethodSelector(row[0].heritability_methods).any { method -> getMethodCapabilities()[method]?.consumes_population_prevalence } ||
                summary_statistics_id in pair_prevalence_summary_statistics_ids
            )
        }
        .collect { row -> normaliseCellValue(row[0].producer_analysis_id).toString() } as Set

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
        validateMethodSelectors(association_methods, heritability_methods, reject, relationship_rows ? true : false)
        def routes = getMethodRoutes(association_methods, heritability_methods)
        def settings = getAnalysisSettings(analysis_meta)
        def is_binary = analysis_meta.trait_type == 'binary'
        validateTraitColumns(is_binary, settings, reject)
        validateMethodConditionedColumns(
            settings,
            routes,
            reject,
            analysis_id in pair_prevalence_analysis_ids || analysis_id in summary_prevalence_analysis_ids,
        )

        if (cohort) {
            def resolved_meta = analysis_meta + [
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
                sample_prevalence: settings.sample_prevalence,
                method_options: method_options_by_analysis[analysis_id],
            ]
            def validated = [
                resolved_meta,
                cohort.genotype_files,
                cells.phenotype,
                cells.quant_covariates,
                cells.cat_covariates,
                method_options_by_analysis[analysis_id].ldak.predictor_extract,
                method_options_by_analysis[analysis_id].ldak.weights,
            ]
            validated_analyses << validated
            analyses_by_id[analysis_id] = validated
            association_methods.each { association_method ->
                def summary_statistics_id = getSummaryStatisticsId(analysis_id, association_method)
                generated_summaries_by_id[summary_statistics_id] = getInternalSummaryMetadata(resolved_meta, association_method)
            }
        }
    }

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
        def source_mode = normaliseCellValue(summary_meta.source_mode)?.toString()
        def source_format = normaliseCellValue(summary_meta.source_format)?.toString()
        def producer_analysis_id = normaliseCellValue(summary_meta.producer_analysis_id)?.toString()
        def producer_association_method = normaliseCellValue(summary_meta.producer_association_method)?.toString()
        def has_external_origin = source || source_mode || source_format
        def has_internal_origin = producer_analysis_id || producer_association_method
        if (has_external_origin && has_internal_origin) {
            reject.call(
                ['source', 'source_mode', 'source_format', 'producer_analysis_id', 'producer_association_method'],
                'external source fields and pipeline-generated producer fields are mutually exclusive',
            )
        }
        if (!has_external_origin && !has_internal_origin) {
            reject.call(
                ['source', 'producer_analysis_id'],
                'no summary-statistics origin is declared; supply a complete external source or producer analysis/method pair',
            )
        }
        if (has_external_origin && (!source || !source_mode || !source_format)) {
            reject.call(['source', 'source_mode', 'source_format'], 'external origin requires all three source fields')
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
        if (has_external_origin && source && source_mode && source_format) {
            if (generated_summaries_by_id.containsKey(summary_statistics_id)) {
                reject.call(
                    'summary_statistics_id',
                    "external result collides with pipeline-generated result '${summary_statistics_id}'; choose a distinct external identity",
                )
            }
            if (source_mode == 'canonical' && source_format != 'nfcore_gwas_canonical_v1') {
                reject.call('source_format', "canonical input must declare 'nfcore_gwas_canonical_v1'")
            }
            if (source_mode == 'raw' && (source_format == 'nfcore_gwas_canonical_v1' || source_format.startsWith('auto'))) {
                reject.call('source_format', "raw input must declare an explicit GWASLab format name and cannot use canonical or automatic detection")
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
                source_mode: source_mode,
                source_format: source_format,
                source_method: summary_meta.source_method,
                source_release: normaliseCellValue(summary_meta.source_release),
                source_name: file(source).name,
                producer_analysis_id: null,
                producer_association_method: null,
                heritability_methods: methods,
                canonical_contract_version: 'nfcore_gwas_canonical_v1',
                transformation: source_mode == 'raw' ? 'gwaslab_harmonised' : 'validated_without_harmonisation',
                harmonization: source_mode == 'raw' ? [tool: 'GWASLab', version: '4.1.9'] : null,
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

        def left_meta = left_analysis ? left_analysis[0] : left_summary
        def right_meta = right_analysis ? right_analysis[0] : right_summary
        if (left_meta && right_meta && left_meta.trait.toString() == right_meta.trait.toString()) {
            reject.call(
                ['left_analysis_id', 'right_analysis_id', 'left_summary_statistics_id', 'right_summary_statistics_id'],
                "declared trait_id '${left_meta.trait}' is equal on both sides; self-pairs are invalid",
            )
        }
        if (analysis_methods && left_analysis && right_analysis && left_analysis[0].cohort.toString() != right_analysis[0].cohort.toString()) {
            reject.call(['left_analysis_id', 'right_analysis_id'], "GCTA bivariate REML requires one cohort, but '${left_analysis_id}' uses '${left_analysis[0].cohort}' and '${right_analysis_id}' uses '${right_analysis[0].cohort}'")
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
                    def matrix_settings = capability.matrix_kind == 'gcta_ldms'
                        ? resolvePairLdmsMatrixSettings(
                            [:],
                            getMethodOptionDefaults().gcta.subMap(['ld_score_region_kb', 'ld_bins', 'ldms_maf_edges']),
                        ) { option, reason -> reject.call('relationship_methods', "default ${option}: ${reason}") }
                        : [:]
                    def pair_meta = common_meta + [
                        cohort: left_analysis_meta.cohort,
                        build: left_analysis_meta.build,
                        ancestry: left_analysis_meta.ancestry,
                        genotype_format: left_analysis_meta.genotype_format,
                        matrix_kind: capability.matrix_kind,
                        matrix_settings: matrix_settings,
                        reml_bivar_prevalence: getGctaBivariatePrevalence(left_analysis_meta, right_analysis_meta),
                    ]
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

    validated_analyses.each { row ->
        def meta = row[0]
        if (!meta.association_methods && !meta.heritability_methods && !(meta.id in referenced_analyses)) {
            def line = line_by_analysis_id[meta.id]
            errors << "  - ${analysis_manifest} row ${line} (analysis_id '${meta.id}'), fields 'association_methods', 'heritability_methods': row selects no unary method and is not referenced by any relationship; populate a method, reference it from a relationship or remove the row"
        }
    }
    declared_summaries.each { meta ->
        if (!meta.heritability_methods && !(meta.summary_statistics_id in referenced_summaries)) {
            def line = line_by_summary_statistics_id[meta.summary_statistics_id]
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
        summary_unary_primary_requests,
        reference_bundles,
    )
    def resolved_pair = resolveRequestNamespace(
        method_options,
        options_document,
        'pair_requests',
        gcta_primary_requests + summary_pair_primary_requests,
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
        analyses: validated_analyses,
        summary_statistics: validated_external_summaries,
        relationships: resolved_relationships,
        unary_requests: unary_requests,
        pair_requests: pair_requests,
    ]
}
