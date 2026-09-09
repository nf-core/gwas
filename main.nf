#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/gwas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/gwas
    Website: https://nf-co.re/gwas
    Slack  : https://nfcore.slack.com/channels/gwas
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GWAS                    } from './workflows/gwas'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_gwas_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_gwas_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.cohort_manifest,
        params.analysis_manifest,
        params.summary_statistics_manifest,
        params.relationship_manifest,
        params.reference_catalog,
        params.method_options,
        params.help,
        params.help_full,
        params.show_hidden,
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_GWAS(
        PIPELINE_INITIALISATION.out.analyses,
        PIPELINE_INITIALISATION.out.summary_statistics,
        PIPELINE_INITIALISATION.out.relationships,
        PIPELINE_INITIALISATION.out.unary_requests,
        PIPELINE_INITIALISATION.out.pair_requests,
        PIPELINE_INITIALISATION.out.meta_requests,
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_GWAS.out.multiqc_report,
    )

    publish:
    genotype_views      = NFCORE_GWAS.out.genotype_views
    gcta_ldms_artifacts = params.save_relatedness_matrices ? NFCORE_GWAS.out.gcta_ldms_artifacts : channel.empty()
    mph_ldms_artifacts  = params.save_relatedness_matrices ? NFCORE_GWAS.out.mph_ldms_artifacts : channel.empty()
    ldms_plan_artifacts = params.save_relatedness_matrices ? NFCORE_GWAS.out.ldms_plan_artifacts : channel.empty()
}

output {
    // Every cohort has a view record explaining its published artifact keys; derived genotype bundles
    // remain work-directory intermediates.
    genotype_views {
        path { cohort_id, _record -> "genotypes/${cohort_id}" }
    }
    gcta_ldms_artifacts {
        path { matrix_meta, _grm_files, _grm_prefixes -> "quality_control/relatedness_matrices/${matrix_meta.key}" }
    }
    mph_ldms_artifacts {
        path { matrix_meta, _grm_files, _grm_prefixes -> "quality_control/relatedness_matrices/${matrix_meta.key}" }
    }
    // The component plan is its own publication family rather than part of either matrix family: two matrix
    // families now consume one plan, so it belongs to neither and is addressed by its own key.
    ldms_plan_artifacts {
        path { plan_meta, _ld_scores, _strata_manifest, _snp_group_files -> "quality_control/ldms_component_plans/${plan_meta.key}" }
    }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_GWAS {
    take:
    analyses // channel: [ val(meta), path(genotype_files), path(phenotype), path(quant_covariates), path(cat_covariates), path(kvik_extract), path(ldak_weights) ]
    summary_statistics // channel: [ val(meta), path(source) ]
    relationships // channel: [ val(meta), path(genotype_files), path(pair_quant_covariates), path(pair_cat_covariates) ]
    unary_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ]
    pair_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ]
    meta_requests // channel: [ val(meta), val(source_summary_statistics_ids) ]

    main:

    //
    // WORKFLOW: Run pipeline
    //
    GWAS(
        analyses,
        summary_statistics,
        relationships,
        unary_requests,
        pair_requests,
        meta_requests,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )

    emit:
    multiqc_report      = GWAS.out.multiqc_report // channel: [ [ path(report) ] ]
    genotype_views      = GWAS.out.genotype_views // channel: [ val(cohort_id), path(genotype_view_record) ], one per cohort
    gcta_ldms_artifacts = GWAS.out.gcta_ldms_artifacts // channel: [ val(matrix_meta), path(grm_files), val(grm_prefixes) ], one per base key
    mph_ldms_artifacts  = GWAS.out.mph_ldms_artifacts // channel: [ val(matrix_meta), path(grm_files), val(grm_prefixes) ], one per base key
    ldms_plan_artifacts = GWAS.out.ldms_plan_artifacts // channel: [ val(plan_meta), path(ld_scores), path(strata_manifest), [ path(snp_group_file), ... ] ], one per plan key
}
