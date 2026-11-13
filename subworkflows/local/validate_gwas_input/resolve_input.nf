include { getSamplesheetPositionalColumns } from './manifest_contracts'
include {
    readMethodOptionsDocument ;
    validateMethodOptions
} from './method_options'
include { resolveAnalyses                 } from './resolve_analyses'
include { resolveCohorts                  } from './resolve_cohorts'

// Validate the lists as a linked contract and construct the canonical tuple consumed by GWAS.
def validateRelationalInput(cohort_rows, analysis_rows, cohort_manifest, analysis_manifest, cohort_schema, analysis_schema, method_options = null) {
    def errors = []
    def options_document = readMethodOptionsDocument(method_options)
    def method_options_by_analysis = validateMethodOptions(method_options, analysis_rows, options_document)

    def cohorts = resolveCohorts(
        cohort_rows,
        getSamplesheetPositionalColumns(cohort_schema),
        cohort_manifest,
        errors,
    )
    def analyses = resolveAnalyses(
        analysis_rows,
        getSamplesheetPositionalColumns(analysis_schema),
        analysis_manifest,
        cohort_manifest,
        cohorts.cohorts_by_id,
        method_options_by_analysis,
        errors,
    )
    if (errors) {
        error("[nf-core/gwas] ERROR: Validation of linked manifests failed!\n\n${errors.join('\n')}\n")
    }

    return [
        analyses: analyses.analyses,
    ]
}
