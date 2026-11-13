/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// MODULE: Local to the pipeline
include { PREPARE_PHENOTYPE_INPUTS        } from '../modules/local/prepare_phenotype_inputs/main'

// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
include { PREPARE_COHORT_GENOTYPES        } from '../subworkflows/local/prepare_cohort_genotypes'
include { ROUTE_GWAS_REPORTING            } from '../subworkflows/local/route_gwas_reporting'

// FUNCTION: Local to the pipeline
include { getMethodTokensWithCapabilities } from '../subworkflows/local/validate_gwas_input/method_registry'

// SUBWORKFLOW: Consisting entirely of nf-core/modules
include { softwareVersionsToYAML          } from '../subworkflows/nf-core/utils_nfcore_pipeline'

// PLUGIN
include { paramsSummaryMap                } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GWAS {
    take:
    ch_analyses // channel: [ val(meta), [ path(genotype_file), ... ], path(phenotype), path(quant_covariates), path(cat_covariates) ]
    multiqc_config // channel: val(multiqc_config)
    multiqc_logo // channel: val(multiqc_logo)
    multiqc_methods_description // channel: val(multiqc_methods_description)
    outdir // channel: val(outdir)

    main:

    // This is the pipeline spine and owns exactly six things: the public `take:` contract above, run-level
    // analysis and method metadata, the genotype consumers and their single per-cohort preparation, phenotype
    // preparation, run-wide version collection and collation, the reporting call, and the public `emit:`
    // block below. Analysis routes attach to the prepared seams; each one is a pipeline-owned subworkflow
    // that receives every configuration value and resource explicitly through its own `take:`, reads neither
    // `params`, `workflow` nor `projectDir`, and owns no shared resource, validation contract or public
    // emission.
    //
    //   ch_analyses
    //        |
    //        v
    //   PREPARE_COHORT_GENOTYPES          PREPARE_PHENOTYPE_INPUTS      shared resources, each built once
    //        |
    //        | cohort_views
    //        v
    //   genotypes/<cohort_id>/<cohort_id>.genotype_view.json
    //
    //   channel.topic('versions') --> softwareVersionsToYAML --> ROUTE_GWAS_REPORTING --> multiqc_report

    //
    // Run-level analysis and method metadata for the report
    //
    // Both are materialised here rather than in the reporting controller because they describe the whole run
    // as validation admitted it, and no single route can see that union.
    def ch_analysis_metadata = ch_analyses
        .map { meta, _genotype_files, _phenotype, _quant_covariates, _cat_covariates -> meta }
        .collect()
    def ch_method_metadata = ch_analyses
        .map { meta, _genotype_files, _phenotype, _quant_covariates, _cat_covariates -> [domain: 'analysis', meta: meta] }
        .collect()

    //
    // The genotype consumers
    //
    // One element per analysis carrying the genotype files it declared. The spine names no representation to
    // decide what is prepared: `PREPARE_COHORT_GENOTYPES` asks the registry which selectors read PLINK 1, and
    // derives that view only for a cohort that has such a consumer. Cohort preparation collapses every request
    // to the distinct cohort before any conversion.
    def ch_genotype_requests = ch_analyses.map { meta, genotype_files, _phenotype, _quant_covariates, _cat_covariates ->
        [meta, genotype_files]
    }

    //
    // SUBWORKFLOW: Prepare each distinct cohort's genotypes once, in the representation it was supplied in
    //
    PREPARE_COHORT_GENOTYPES(ch_genotype_requests)

    //
    // MODULE: Prepare each analysis's phenotype and covariates once
    //
    // Which selectors read a covariate file through an interface that treats a missing cell as a value is a
    // registry capability, so the answer is resolved once here rather than re-derived per row. It belongs in
    // the preparation meta, and therefore inside the preparation cache boundary, precisely because it decides
    // whether preparation succeeds at all.
    def complete_covariate_methods = getMethodTokensWithCapabilities([requires_complete_covariates: true])
    PREPARE_PHENOTYPE_INPUTS(
        ch_analyses.map { meta, _genotype_files, phenotype, quant_covariates, cat_covariates ->
            def preparation_meta = [
                id: meta.id,
                phenotype_column: meta.phenotype_column,
                is_binary: meta.is_binary,
                case_value: meta.case_value,
                control_value: meta.control_value,
                requires_complete_covariates: (meta.association_methods + meta.heritability_methods).any { method -> method in complete_covariate_methods },
            ]
            [preparation_meta, phenotype, quant_covariates, cat_covariates]
        }
    )

    //
    // Collate and save software versions
    //
    def ch_topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def ch_topic_versions_string = ch_topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_topic_versions.versions_file)
        .mix(ch_topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_' + 'gwas_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    //
    // SUBWORKFLOW: Render the run report from the analysis plan, workflow summary, methods description and versions
    //
    // The reporting controller owns MultiQC assembly but reads no parent scope: the run parameter summary is
    // evaluated here and every pipeline-default asset is resolved here, then passed in explicitly.
    def summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ROUTE_GWAS_REPORTING(
        ch_collated_versions,
        ch_analysis_metadata,
        ch_method_metadata,
        summary_params,
        multiqc_config,
        multiqc_logo,
        multiqc_methods_description,
        file("${projectDir}/assets/multiqc_analysis_plan.yml", checkIfExists: true),
        file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
        file("${projectDir}/assets/nf-core-gwas_logo_light.png", checkIfExists: true),
        file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true),
    )

    // One JSON provenance record per cohort. It is the document that explains every published artifact key a
    // run produced — which bytes the cohort's identity was taken from, which representation was actually on
    // disk, and, where a PLINK 1 view was derived, the projection policy and what that policy discarded — so
    // it is written unconditionally rather than behind a save control. `cohort_views` is its only consumer,
    // so the late emission of the joined record it is built from affects nothing else.
    def ch_genotype_view_records = PREPARE_COHORT_GENOTYPES.out.cohort_views
        .collectFile { cohort_meta, view ->
            ["${cohort_meta.id}.genotype_view.json", groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(view)) + '\n']
        }
        .map { record -> [record.name - '.genotype_view.json', record] }

    emit:
    multiqc_report = ROUTE_GWAS_REPORTING.out.report.toList() // channel: [ [ path(report) ] ]
    genotype_views = ch_genotype_view_records // channel: [ val(cohort_id), path(genotype_view_record) ], one per cohort
}
