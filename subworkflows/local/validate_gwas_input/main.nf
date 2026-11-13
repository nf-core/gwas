// Validate linked cohort and analysis manifests as one pipeline input contract and emit canonical analyses.
// This workflow performs no tasks and therefore emits no versions.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// PLUGIN
include { samplesheetToList                } from 'plugin/nf-schema'

// LOCAL FUNCTIONS
include { validateUniqueSamplesheetHeaders } from './manifest_contracts'
include { validateRelationalInput          } from './resolve_input'

workflow VALIDATE_GWAS_INPUT {
    take:
    cohort_manifest // channel: val(cohort_manifest)
    analysis_manifest // channel: val(analysis_manifest)
    method_options // channel: val(method_options)

    main:
    if (analysis_manifest && !cohort_manifest) {
        error("[nf-core/gwas] ERROR: the analysis manifest --input '${analysis_manifest ?: ''}' was supplied but --cohort_manifest is missing; linked manifests require both")
    }
    if (cohort_manifest && !analysis_manifest) {
        error("[nf-core/gwas] ERROR: --cohort_manifest '${cohort_manifest}' was supplied but the analysis manifest --input is missing; linked manifests require both")
    }
    if (!cohort_manifest && !analysis_manifest) {
        error("[nf-core/gwas] ERROR: supply the analysis manifest as --input and the linked cohort manifest as --cohort_manifest")
    }

    def cohort_schema = "${projectDir}/assets/schema_cohort_manifest.json"
    def analysis_schema = "${projectDir}/assets/schema_analysis_manifest.json"
    validateUniqueSamplesheetHeaders(cohort_manifest, 'Cohort manifest')
    validateUniqueSamplesheetHeaders(analysis_manifest, 'Analysis manifest')

    def validated = validateRelationalInput(
        samplesheetToList(cohort_manifest, cohort_schema),
        samplesheetToList(analysis_manifest, analysis_schema),
        cohort_manifest,
        analysis_manifest,
        cohort_schema,
        analysis_schema,
        method_options,
    )
    def ch_analyses = channel.fromList(validated.analyses)

    emit:
    analyses = ch_analyses // channel: [ val(meta), [ path(genotype_file), ... ], path(phenotype), path(quant_covariates), path(cat_covariates) ]
}
