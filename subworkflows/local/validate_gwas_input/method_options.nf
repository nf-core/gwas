include { tokenizeMethodSelector          } from './manifest_contracts'

// The accepted method families, their accepted option names and every resolved default are this one map.
// Family acceptance, per-family key validation and the resolved shape are all derived from it, so adding a
// family is adding an entry here rather than editing three parallel lists. Families are keyed alphabetically.
// No method is registered yet, so the document accepts no family: every analysis resolves to the empty map.
def getMethodOptionDefaults() {
    return [:]
}

// Operational tuning that a researcher might reasonably try to place in the scientific option document. Naming
// it per family lets the diagnostic say "run/profile configuration" instead of "unknown option".
def getMethodOptionOperationalKeys() {
    return [:]
}

def validateMethodOptionKeys(analysis_id, family, options, accepted, operational, fail) {
    def unknown = options.keySet().findAll { option -> !(option in accepted) }
    if (unknown) {
        def option = unknown.first()
        def reason = option in operational
            ? 'operational tuning must be supplied through run/profile configuration'
            : "unknown option; accepted ${family.toUpperCase()} options are ${accepted.join(', ')}"
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
    def namespaces = ['analyses']
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
        // so a family's accepted names and its resolved shape come from one place.
        def unknown_families = families.keySet().findAll { family -> !defaults.containsKey(family) }
        if (unknown_families) {
            def accepted = defaults ? defaults.keySet().join(', ') : 'none, no method is registered in this release'
            fail.call(analysis_id, unknown_families.first().toString(), "unknown method family; accepted families are ${accepted}")
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

        resolved[analysis_id] = defaults.keySet().collectEntries { family -> [(family): supplied[family]] }
    }

    return resolved
}
