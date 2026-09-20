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
    catch (exception: Exception) {
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
            def optional = []
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
            // LDSC's `--ref-ld-chr` and `--w-ld-chr` take a stem, not a directory: LDSC appends
            // `<chr>.l2.ldscore.gz` to it, and `<chr>.l2.M_5_50` for the reference set. The catalog records
            // the text a person would type — `/refs/eur_w_ld_chr/` or `/refs/baseline/baselineLD.` — and the
            // split at the last `/` recovers the directory the pipeline stages and the file-name prefix the
            // LDSC modules hand back to the tool. The declared text is authoritative; nothing is inferred
            // from the directory listing.
            def stem_fields = family == 'ldsc' ? ['reference_ld_scores', 'regression_weights'] : []
            def resources = [:]
            def stem_prefixes = [:]
            def role_names = [:]
            role_fields.each { field ->
                def declared = definition[field].toString()
                if (field in stem_fields) {
                    def separator = declared.lastIndexOf('/')
                    def directory_text = separator < 0 ? '.' : separator == 0 ? '/' : declared.substring(0, separator)
                    def prefix = separator < 0 ? declared : declared.substring(separator + 1)
                    def directory = file(directory_text)
                    if (!directory.exists()) {
                        fail.call(bundle_id, field, "stem '${declared}' names directory '${directory_text}', which does not exist")
                    }
                    if (!directory.toFile().isDirectory()) {
                        fail.call(bundle_id, field, "stem '${declared}' names directory '${directory_text}', which must be a directory")
                    }
                    def expected = field == 'reference_ld_scores'
                        ? ["${prefix}1.l2.ldscore.gz", "${prefix}1.l2.M_5_50"]
                        : ["${prefix}1.l2.ldscore.gz"]
                    expected.each { name ->
                        def candidate = directory.resolve(name.toString())
                        if (!candidate.exists()) {
                            fail.call(bundle_id, field, "stem '${declared}' does not resolve chromosome 1; expected '${candidate}'")
                        }
                    }
                    resources[field] = directory
                    stem_prefixes[field] = prefix
                    role_names[field] = "${directory.name}/${prefix}".toString()
                }
                else {
                    def resource = file(declared)
                    if (!resource.exists()) {
                        fail.call(bundle_id, field, "resource path '${declared}' does not exist")
                    }
                    if (!resource.toFile().isFile()) {
                        fail.call(bundle_id, field, "resource path '${declared}' must be a file")
                    }
                    resources[field] = resource
                    role_names[field] = resource.name
                }
            }
            resolved[bundle_id] = [
                id: bundle_id,
                family: family,
                genome_build: definition.genome_build,
                ancestry: definition.ancestry,
                variant_id_system: definition.variant_id_system,
                model: family == 'ldak' ? definition.model : null,
                role_names: role_names,
                resources: resources,
                stem_prefixes: stem_prefixes,
            ]
        }
    }
    return resolved
}
