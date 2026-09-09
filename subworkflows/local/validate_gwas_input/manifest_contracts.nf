include { getAssociationMethodTokens ; getHeritabilityMethodTokens ; getMethodCapabilities ; getMethodTokensWithCapabilities } from './method_registry'

def getGenotypeGroups() {
    return [
        plink2: ['pgen', 'psam', 'pvar'],
        plink1: ['bed', 'bim', 'fam'],
    ]
}

// Derive nf-schema positional fields from each schema so adding a property cannot silently shift the tuple.
def getSamplesheetPositionalColumns(schema) {
    def properties = new groovy.json.JsonSlurper().parseText(file(schema).text).items.properties
    return properties.findAll { _column, definition -> !definition.containsKey('meta') }.keySet().toList()
}

// nf-schema applies the schema contract but cannot distinguish repeated names once it has created a map.
// Read the header as CSV solely to preserve that one piece of ingress information for validation.
def validateUniqueSamplesheetHeaders(samplesheet, role) {
    def header_line = file(samplesheet).readLines().find { line -> line.trim() }
    def observed = header_line
        ? header_line.split(/,(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)/, -1).collect { name -> name.trim().replaceAll(/^\"|\"$/, '').replace('""', '"') }
        : []
    def repeated = observed.countBy { column -> column }.findAll { _column, count -> count > 1 }.keySet().toList()

    if (repeated) {
        error("[nf-core/gwas] ERROR: ${role} '${samplesheet}' header row 1 repeats column${repeated.size() > 1 ? 's' : ''} ${repeated.collect { name -> "'${name}'" }.join(', ')}.\n")
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
    def ldak_association = getMethodTokensWithCapabilities([domain: 'association', option_family: 'ldak'])
    def sparse_grm_association = getMethodTokensWithCapabilities([domain: 'association', input_backend: 'sparse_grm'])
    def ldms_heritability = getMethodTokensWithCapabilities([domain: 'heritability', input_backend: 'ldms_grm_family'])
    return [
        runs_heritability: selected.any { details -> details.domain == 'heritability' },
        consumes_population_prevalence: selected.any { details -> details.prevalence.population != 'not_consumed' },
        runs_ldak_kvik: association_methods.any { method -> method in ldak_association },
        runs_ldak_heritability: heritability_methods.any { method -> capabilities[method] && capabilities[method].option_family == 'ldak' },
        requires_population_prevalence: heritability_methods.findAll { method -> capabilities[method] && capabilities[method].prevalence.population == 'required' },
        runs_gcta: selected.any { details -> details.option_family == 'gcta' },
        runs_gcta_fastgwa: association_methods.any { method -> method in sparse_grm_association },
        runs_greml_ldms: heritability_methods.any { method -> method in ldms_heritability },
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

// Enforce each estimator's declared trait support at manifest ingress. The PCGC modes require a binary
// phenotype; MPH has no liability-scale conversion. The diagnostic names the trait type and capable methods.
def validateHeritabilityTraitSupport(is_binary, heritability_methods, reject) {
    def declared = is_binary ? 'binary' : 'quantitative'
    def capabilities = getMethodCapabilities()
    def unsupported = heritability_methods
        .findAll { method ->
            capabilities[method] && !capabilities[method].trait_support[declared]
        }
        .unique()
    if (!unsupported) {
        return null
    }
    def capable = getHeritabilityMethodTokens().findAll { method -> capabilities[method].trait_support[declared] }
    reject.call(
        'heritability_methods',
        "method(s) ${unsupported.join(', ')} support ${is_binary ? 'quantitative' : 'binary'} traits only, but this row declares trait_type '${declared}'; ${declared}-capable heritability methods are ${capable.join(', ')}",
    )
}

def validateGenotypeGroup(cells, reject) {
    def groups = getGenotypeGroups()
    def populated_groups = groups.findAll { _name, columns -> columns.any { column -> cells[column] } }
    if (!populated_groups) {
        reject.call(groups.values().flatten(), 'no genotype group is populated, supply exactly one of pgen/psam/pvar or bed/bim/fam')
    }
    else if (populated_groups.size() > 1) {
        populated_groups
            .keySet()
            .toList()
            .tail()
            .each { name ->
                def populated_column = groups[name].find { column -> cells[column] }
                reject.call(populated_column, 'a second genotype group is populated on this row, supply exactly one of pgen/psam/pvar or bed/bim/fam')
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

// Every PLINK consumer in the pipeline addresses a fileset by one prefix rather than by three paths:
// `plink2 --pfile`, `plink2 --bfile`, `regenie --bed`, `gcta --mbfile`/`--mpfile` and every LDAK executable
// all take the shared basename stem and append the member extensions themselves. The manifest, by contrast,
// declares each member as an independent `file-path` property, so a trio whose members do not share one stem
// passes the schema and is only discovered by whichever consumer runs first, each with a different message
// and none naming the cohort. Rejecting it here is ingress validation of a user-supplied row, and it is the
// one place that can name the cohort and the differing columns.
def validateGenotypeBasenameStem(genotype_format, cells, reject) {
    if (!genotype_format) {
        return null
    }
    def columns = getGenotypeGroups()[genotype_format]
    // A single-member group has no stem to agree on, and an incomplete group has already been rejected by
    // `validateGenotypeGroup`, so neither is reported a second time here.
    if (columns.size() < 2 || columns.any { column -> !cells[column] }) {
        return null
    }
    def stems = columns.collectEntries { column -> [(column): file(cells[column].toString()).baseName] }
    if (stems.values().toList().unique().size() > 1) {
        reject.call(
            columns,
            "genotype group '${genotype_format}' members do not share one basename stem (${columns.collect { column -> "${column}='${stems[column]}'" }.join(', ')}); every PLINK, GCTA and REGENIE consumer addresses the fileset by a single prefix, so the three members must be named <stem>.${columns.join(', <stem>.')}",
        )
    }
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
    if (settings.population_prevalence == null && routes.requires_population_prevalence) {
        reject.call('population_prevalence', "method(s) ${routes.requires_population_prevalence.join(', ')} always estimate on the liability scale and require a population prevalence")
    }
}
