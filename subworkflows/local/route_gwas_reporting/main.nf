// Assemble the nf-core/gwas run report: the analysis plan, the workflow summary, the route-aware methods
// description and the run-wide collated software versions, then render them once with MultiQC.
// This is pipeline reporting policy, not an nf-core/modules submission candidate.
// Every configuration value and asset arrives through `take:`; nothing here reads params or projectDir.

// MODULE: Installed directly from nf-core/modules
include { MULTIQC                } from '../../../modules/nf-core/multiqc/main'

// FUNCTION: Local to the pipeline
include { analysisPlanJson       } from '../utils_nfcore_gwas_pipeline'
include { methodsDescriptionText } from '../utils_nfcore_gwas_pipeline'

// FUNCTION: Installed directly from nf-core/subworkflows
include { paramsSummaryMultiqc   } from '../../nf-core/utils_nfcore_pipeline'

workflow ROUTE_GWAS_REPORTING {
    take:
    ch_collated_versions // channel: path(collated_versions)
    ch_analysis_metadata // channel: val([ val(meta), ... ]), every validated analysis unit, collected once
    ch_method_metadata // channel: val([ val([domain, meta]), ... ]), every method-selecting record, collected once
    summary_params // channel: val(summary_params), the evaluated run parameter summary grouped by schema section
    multiqc_config // channel: val(multiqc_config), user-supplied override or null
    multiqc_logo // channel: val(multiqc_logo), user-supplied override or null
    multiqc_methods_description // channel: val(multiqc_methods_description), user-supplied override or null
    multiqc_analysis_plan // channel: val(multiqc_analysis_plan), resolved analysis-plan section template
    default_multiqc_config // channel: val(default_multiqc_config), resolved pipeline default
    default_multiqc_logo // channel: val(default_multiqc_logo), resolved pipeline default
    default_methods_description // channel: val(default_methods_description), resolved pipeline default

    main:

    // Every section is rendered once from an already-collected metadata view, so each contributes exactly
    // one file to the single MultiQC task. Route controllers contribute nothing else to this channel.
    def ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    def ch_analysis_plan = ch_analysis_metadata.map { analysis_metadata -> analysisPlanJson(multiqc_analysis_plan, analysis_metadata) }
    ch_multiqc_files = ch_multiqc_files.mix(ch_analysis_plan.collectFile(name: 'analysis_plan_mqc.json', sort: true))

    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    def multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : default_methods_description
    def ch_methods_description = ch_method_metadata.map { method_metadata ->
        methodsDescriptionText(multiqc_custom_methods_description, getSelectedMethods(method_metadata))
    }
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    //
    // MODULE: MultiQC
    //
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'gwas'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : default_multiqc_config,
                multiqc_logo
                    ? file(multiqc_logo, checkIfExists: true)
                    : default_multiqc_logo,
                [],
                [],
            ]
        }
    )

    emit:
    report = MULTIQC.out.report.map { _meta, report -> [report] } // channel: [ path(multiqc_report) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Reduce the collected method-selecting records to the distinct method tokens the run actually requested,
// grouped by the domain the methods description cites. Summary-scale unary requests name a single method
// and are cited alongside the individual-level heritability methods.
def getSelectedMethods(method_metadata) {
    def analysis_metadata = method_metadata.findAll { record -> record.domain == 'analysis' }.collect { record -> record.meta }
    def summary_unary_metadata = method_metadata.findAll { record -> record.domain == 'summary_unary' }.collect { record -> record.meta }
    def relationship_metadata = method_metadata.findAll { record -> record.domain == 'pairwise' }.collect { record -> record.meta }
    def summary_set_metadata = method_metadata.findAll { record -> record.domain == 'summary_set' }.collect { record -> record.meta }
    return [
        association: analysis_metadata.collectMany { meta -> meta.association_methods }.unique().sort(),
        heritability: (analysis_metadata.collectMany { meta -> meta.heritability_methods } + summary_unary_metadata.collect { meta -> meta.method }).unique().sort(),
        pairwise: relationship_metadata.collectMany { meta -> meta.relationship_methods }.unique().sort(),
        summary_set: summary_set_metadata.collectMany { meta -> meta.meta_analysis_models }.unique().sort(),
    ]
}
