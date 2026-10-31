include { getGenotypeGroups ; normaliseCellValue ; validateGenotypeBasenameStem ; validateGenotypeGroup } from './manifest_contracts'

def resolveCohorts(cohort_rows, cohort_columns, cohort_manifest, errors) {
    def cohorts_by_id = [:]
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
        validateGenotypeBasenameStem(genotype_format, cells, reject)
        def genotype_files = genotype_format
            ? getGenotypeGroups()[genotype_format].collect { column -> cells[column] }
            : []

        // The declared view identity joins the definition rather than sitting beside it, so two rows that
        // share a `cohort_id` but declare different identities conflict by name here instead of silently
        // taking whichever row was read first. It is normalised to '' when absent so an absent and a blank
        // cell are the same definition.
        def genotype_view_id = normaliseCellValue(cells.genotype_view_id)?.toString()
        def definition = [
            genome_build: cohort_meta.build,
            ancestry: cohort_meta.ancestry,
            genotype_view_id: genotype_view_id ?: '',
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
                genotype_view_id: genotype_view_id,
                definition: definition,
            ]
        }
    }
    return [cohorts_by_id: cohorts_by_id]
}
