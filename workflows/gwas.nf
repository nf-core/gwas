include { ROUTE_META_ANALYSIS } from '../subworkflows/local/route_meta_analysis'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// MODULE: Local to the pipeline
include { PREPARE_PHENOTYPE_INPUTS           } from '../modules/local/prepare_phenotype_inputs/main'

// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
include { PREPARE_COHORT_GENOTYPES           } from '../subworkflows/local/prepare_cohort_genotypes'
include { PREPARE_RELATEDNESS_MATRICES       } from '../subworkflows/local/prepare_relatedness_matrices'
include { PREPARE_RELATIONSHIP_TRAITS        } from '../subworkflows/local/prepare_relationship_traits'
include { ROUTE_ASSOCIATION_ANALYSES         } from '../subworkflows/local/route_association_analyses'
include { ROUTE_GRM_HERITABILITY             } from '../subworkflows/local/route_grm_heritability'
include { ROUTE_LDAK_DIRECT_HERITABILITY     } from '../subworkflows/local/route_ldak_direct_heritability'
include { ROUTE_MPH_HERITABILITY             } from '../subworkflows/local/route_mph_heritability'
include { ROUTE_GCTA_BIVARIATE_RELATIONSHIPS } from '../subworkflows/local/route_gcta_bivariate_relationships'
include { ROUTE_MPH_BIVARIATE_RELATIONSHIPS  } from '../subworkflows/local/route_mph_bivariate_relationships'
include { ROUTE_CANONICAL_SUMMARY_STATISTICS } from '../subworkflows/local/route_canonical_summary_statistics'
include { ROUTE_LDAK_SUMMARY_ANALYSES        } from '../subworkflows/local/route_ldak_summary_analyses'
include { ROUTE_LDSC_SUMMARY_ANALYSES        } from '../subworkflows/local/route_ldsc_summary_analyses'
include { ROUTE_GWAS_REPORTING               } from '../subworkflows/local/route_gwas_reporting'
include { getGwaslabReferences               } from '../subworkflows/local/utils_nfcore_gwas_pipeline'

// FUNCTION: Local to the pipeline
include { getMethodTokensWithCapabilities    } from '../subworkflows/local/validate_gwas_input/method_registry'

// SUBWORKFLOW: Consisting entirely of nf-core/modules
include { softwareVersionsToYAML             } from '../subworkflows/nf-core/utils_nfcore_pipeline'

// PLUGIN
include { paramsSummaryMap                   } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GWAS {
    take:
    ch_analyses // channel: [ val(meta), [ path(genotype_file), ... ], path(phenotype), path(quant_covariates), path(cat_covariates), path(kvik_extract), path(ldak_weights) ]
    ch_external_summary_statistics // channel: [ val(meta), path(source) ]
    ch_relationships // channel: [ val(meta), [ path(genotype_file), ... ], path(pair_quant_covariates), path(pair_cat_covariates) ]
    ch_unary_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), val(reference_prefix), path(regression_weights), val(weights_prefix), path(tagging_file) ]
    ch_pair_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), val(reference_prefix), path(regression_weights), val(weights_prefix), path(tagging_file) ]
    ch_meta_requests // channel: [ val(meta), val(source_summary_statistics_ids) ]
    multiqc_config // channel: val(multiqc_config)
    multiqc_logo // channel: val(multiqc_logo)
    multiqc_methods_description // channel: val(multiqc_methods_description)
    outdir // channel: val(outdir)

    main:

    // This is the pipeline spine and owns exactly ten things: the public `take:` contract above, run-level
    // analysis and method metadata, the union of genotype consumers and their single preparation, the union
    // of relatedness-matrix consumers and their single preparation, phenotype preparation plus the
    // tool-neutral per-analysis seams derived from it, the route-controller calls and the dependencies
    // between their semantic results, the fan-out of GWASLab-standard summaries to the summary-scale routes,
    // run-wide version collection and collation, the reporting call, and the public `emit:` block below.
    //
    // Every route controller is a pipeline-owned subworkflow that receives all configuration values and
    // resources explicitly through its own `take:`. None of them reads `params`, `workflow` or `projectDir`,
    // and none of them owns a shared resource, the validation contract, or a public emission.
    //
    //   ch_analyses / ch_relationships
    //          |
    //          v
    //   PREPARE_COHORT_GENOTYPES ---> PREPARE_RELATEDNESS_MATRICES     PREPARE_PHENOTYPE_INPUTS
    //          |          |                     |                              |
    //          |          |                     |                    PREPARE_RELATIONSHIP_TRAITS
    //          +----------+---------------------+------------------------------+   shared resources,
    //          |          |                     |                              |   each built once
    //          v          v                     v                              v
    //   ROUTE_ASSOCIATION_ANALYSES     ROUTE_GRM_HERITABILITY     ROUTE_GCTA_BIVARIATE_RELATIONSHIPS
    //          |          |                                       ROUTE_MPH_BIVARIATE_RELATIONSHIPS
    //          |   ROUTE_MPH_HERITABILITY           (native MPH matrices built on the spine beside the others)
    //          |   ROUTE_LDAK_DIRECT_HERITABILITY   (direct genotypes; requests no relatedness matrix)
    //          |
    //          | association_results                        ch_external_summary_statistics
    //          v                                                         |
    //   ROUTE_CANONICAL_SUMMARY_STATISTICS <---------------------------- +
    //          |
    //          | summary_statistics (GWASLab convergence point)
    //          +--> ROUTE_LDAK_SUMMARY_ANALYSES
    //          +--> ROUTE_LDSC_SUMMARY_ANALYSES
    //          +--> ROUTE_META_ANALYSIS combines the canonical study summaries
    //
    //   channel.topic('versions') --> softwareVersionsToYAML --> ROUTE_GWAS_REPORTING --> multiqc_report

    //
    // Run-level analysis and method metadata for the report
    //
    // Both are materialised here rather than in the reporting controller because they describe the whole run
    // as validation admitted it, across all four request domains, and no single route can see that union.
    def ch_analysis_metadata = ch_analyses
        .map { meta, _genotype_files, _phenotype, _quant_covariates, _cat_covariates, _kvik_extract, _ldak_weights -> meta }
        .collect()
    def ch_method_metadata = ch_analyses
        .map { meta, _genotype_files, _phenotype, _quant_covariates, _cat_covariates, _kvik_extract, _ldak_weights -> [domain: 'analysis', meta: meta] }
        .mix(
            ch_relationships.map { meta, _genotype_files, _pair_quant_covariates, _pair_cat_covariates -> [domain: 'pairwise', meta: meta] }
        )
        .mix(
            ch_unary_requests.map { meta, _hapmap3_snplist, _reference_ld_scores, _reference_prefix, _regression_weights, _weights_prefix, _tagging_file -> [domain: 'summary_unary', meta: meta] }
        )
        .mix(
            ch_pair_requests.map { meta, _hapmap3_snplist, _reference_ld_scores, _reference_prefix, _regression_weights, _weights_prefix, _tagging_file -> [domain: 'pairwise', meta: meta] }
        )
        .mix(
            ch_meta_requests.map { meta, _source_ids -> [domain: 'summary_set', meta: meta] }
        )
        .collect()

    //
    // Union of the genotype consumers across every request domain
    //
    // One element per analysis unit carrying the genotype files it declared. Every pairwise request is also a
    // potential consumer of the cohort's PLINK 1 view through its matrix kind's declared genotype bundle, so
    // it enters this request stream without inheriting either endpoint's unary method settings. The spine
    // names no matrix kind and no representation to decide that: `PREPARE_COHORT_GENOTYPES` asks the registry
    // which selectors read PLINK 1, and derives the view only for a cohort that has such a consumer. Cohort
    // preparation collapses every request to the distinct cohort before any conversion.
    def ch_genotype_requests = ch_analyses.map { meta, genotype_files, _phenotype, _quant_covariates, _cat_covariates, _kvik_extract, _ldak_weights ->
        [meta, genotype_files]
    }
    ch_genotype_requests = ch_genotype_requests.mix(
        ch_relationships.map { meta, genotype_files, _pair_quant_covariates, _pair_cat_covariates -> [meta, genotype_files] }
    )

    //
    // Union of the relatedness-matrix consumers across every request domain
    //
    // Matrix preparation additionally receives the optional LDAK weights Path. It derives identity from the
    // bytes before request deduplication and keeps the Path outside matrix metadata and the published key.
    def ch_relatedness_analyses = ch_analyses.map { meta, genotype_files, _phenotype, _quant_covariates, _cat_covariates, _kvik_extract, ldak_weights ->
        [meta, genotype_files, ldak_weights ?: []]
    }
    ch_relatedness_analyses = ch_relatedness_analyses.mix(
        ch_relationships.map { meta, genotype_files, _pair_quant_covariates, _pair_cat_covariates -> [meta, genotype_files, []] }
    )

    //
    // SUBWORKFLOW: Prepare each distinct cohort's genotypes once, in the representation it was supplied in
    //
    PREPARE_COHORT_GENOTYPES(ch_genotype_requests)

    //
    // SUBWORKFLOW: Build each distinct relatedness matrix once and fan it out per analysis unit
    //
    PREPARE_RELATEDNESS_MATRICES(
        ch_relatedness_analyses,
        PREPARE_COHORT_GENOTYPES.out.cohort_native_genotypes,
        PREPARE_COHORT_GENOTYPES.out.plink1_genotypes,
        PREPARE_COHORT_GENOTYPES.out.cohort_native_view_keys,
        PREPARE_COHORT_GENOTYPES.out.cohort_plink1_view_keys,
        params.gcta_grm_parts,
    )

    //
    // MODULE: Prepare each analysis unit's phenotype and covariates in the shared tool-compatible layout
    //
    // The preparation task receives only the fields its template and configured prefix consume. The complete
    // focal analysis map is restored after the task, so downstream method, request and display metadata cannot
    // alter the preparation cache boundary while every consumer still receives its original analysis identity.
    def ch_analysis_meta_by_id = ch_analyses.map { meta, _genotype_files, _phenotype, _quant_covariates, _cat_covariates, _kvik_extract, _ldak_weights ->
        [meta.id, meta]
    }

    // Which selectors read a covariate file through an interface that treats a missing cell as a value is a
    // registry capability, so the answer is resolved once here rather than re-derived per row. It belongs in
    // the preparation meta, and therefore inside the preparation cache boundary, precisely because it decides
    // whether preparation succeeds at all. LDAK 6.1
    // (pinned genomedk build) fits a stale read-buffer value for a missing --covar cell at exit 0.
    // Retire when the LDAK analysis-row modules move to 6.3, which mean-imputes as documented.
    def complete_covariate_methods = getMethodTokensWithCapabilities([requires_complete_covariates: true])
    PREPARE_PHENOTYPE_INPUTS(
        ch_analyses.map { meta, _genotype_files, phenotype, quant_covariates, cat_covariates, _kvik_extract, _ldak_weights ->
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

    def ch_prepared_phenotype = PREPARE_PHENOTYPE_INPUTS.out.phenotype
        .map { preparation_meta, phenotype -> [preparation_meta.id, phenotype] }
        .join(ch_analysis_meta_by_id, failOnMismatch: true, failOnDuplicate: true)
        .map { _analysis_id, phenotype, meta -> [meta, phenotype] }
    def ch_prepared_phenotype_headerless = PREPARE_PHENOTYPE_INPUTS.out.phenotype_headerless
        .map { preparation_meta, phenotype -> [preparation_meta.id, phenotype] }
        .join(ch_analysis_meta_by_id, failOnMismatch: true, failOnDuplicate: true)
        .map { _analysis_id, phenotype, meta -> [meta, phenotype] }
    def ch_prepared_quant_covariates_headerless = PREPARE_PHENOTYPE_INPUTS.out.quant_covariates_headerless
        .map { preparation_meta, covariates -> [preparation_meta.id, covariates] }
        .join(ch_analysis_meta_by_id, failOnDuplicate: true)
        .map { _analysis_id, covariates, meta -> [meta, covariates] }
    def ch_prepared_cat_covariates_headerless = PREPARE_PHENOTYPE_INPUTS.out.cat_covariates_headerless
        .map { preparation_meta, covariates -> [preparation_meta.id, covariates] }
        .join(ch_analysis_meta_by_id, failOnDuplicate: true)
        .map { _analysis_id, covariates, meta -> [meta, covariates] }
    def ch_prepared_covariates = PREPARE_PHENOTYPE_INPUTS.out.covariates
        .map { preparation_meta, covariates -> [preparation_meta.id, covariates] }
        .join(ch_analysis_meta_by_id, failOnDuplicate: true)
        .map { _analysis_id, covariates, meta -> [meta, covariates] }
    def ch_prepared_adjustment_covariates = PREPARE_PHENOTYPE_INPUTS.out.adjustment_covariates
        .map { preparation_meta, covariates -> [preparation_meta.id, covariates] }
        .join(ch_analysis_meta_by_id, failOnDuplicate: true)
        .map { _analysis_id, covariates, meta -> [meta, covariates] }

    // GCTA and LDAK reject a header row. LDAK-KVIK, fastGWA and every individual-level GRM heritability
    // estimator therefore consume the headerless phenotype and covariate serialisations. This one prepared
    // stream is built here because it has consumers in more than one route, and is passed to the association
    // and heritability controllers explicitly. Optional covariates are represented by [], which stages nothing.
    def ch_gcta_phenotypes = ch_prepared_phenotype_headerless
        .join(ch_prepared_quant_covariates_headerless, remainder: true)
        .join(ch_prepared_cat_covariates_headerless, remainder: true)
        .map { meta, phenotype, quant_covariates, cat_covariates ->
            [meta, phenotype, quant_covariates ?: [], cat_covariates ?: []]
        }

    // MPH's covariate interface is name-keyed: it takes a comma-separated list of column names and reports the
    // fitted effects under those names, so its serializer consumes the headered prepared tables rather than
    // the headerless serialisations GCTA and LDAK read. This seam is built here rather than in the route
    // because deriving per-analysis streams from preparation is spine work, and one absent covariate table is
    // [] so it stages nothing.
    def ch_prepared_covariate_tables = ch_analysis_meta_by_id
        .join(PREPARE_PHENOTYPE_INPUTS.out.quant_covariates.map { preparation_meta, covariates -> [preparation_meta.id, covariates] }, remainder: true)
        .join(PREPARE_PHENOTYPE_INPUTS.out.cat_covariates.map { preparation_meta, covariates -> [preparation_meta.id, covariates] }, remainder: true)
        .map { _analysis_id, meta, quant_covariates, cat_covariates -> [meta, quant_covariates ?: [], cat_covariates ?: []] }

    //
    // SUBWORKFLOW: Prepare each relationship's ordered two-trait table and pair covariates once
    //
    // A relationship is prepared once however many individual-level pair methods select it: the ordered union
    // table and the normalised pair covariates are the scientific pair, and every method adapter serialises
    // that one artifact rather than resolving the endpoints again. Both pair controllers below consume it.
    PREPARE_RELATIONSHIP_TRAITS(ch_relationships, PREPARE_PHENOTYPE_INPUTS.out.phenotype_headerless)

    //
    // SUBWORKFLOW: Pipeline route for REGENIE, LDAK-KVIK and GCTA fastGWA associations
    //
    // Cohort genotype preparation, relatedness-matrix construction and phenotype preparation stay above on
    // the spine so each shared resource is built once and fanned out to every consumer across every domain.
    // The controller owns association-method selection, the adaptation of those prepared streams into each
    // family's native call shape, the two prediction-reusing routes, and the fan-in of three native result
    // contracts onto one raw-association stream naming the producing method.
    //
    // The declared optional LDAK predictor list is narrowed out of the validated relational row here, because
    // reading that row is spine knowledge; which analyses want it, and what its content identity contributes
    // to the Step 1 reuse key, is the controller's. An absent file is [] and stages nothing.
    ROUTE_ASSOCIATION_ANALYSES(
        PREPARE_COHORT_GENOTYPES.out.native_genotypes,
        PREPARE_COHORT_GENOTYPES.out.plink1_genotypes,
        PREPARE_COHORT_GENOTYPES.out.cohort_plink1_view_keys,
        ch_prepared_phenotype,
        ch_prepared_covariates,
        ch_gcta_phenotypes,
        PREPARE_RELATEDNESS_MATRICES.out.gcta_sparse,
        ch_analyses.map { meta, _genotype_files, _phenotype, _quant_covariates, _cat_covariates, kvik_extract, _ldak_weights ->
            [meta, kvik_extract ?: []]
        },
        params.regenie_step2_bsize,
        params.regenie_step1_mode,
        params.regenie_step1_jobs,
    )

    //
    // SUBWORKFLOW: Pipeline route for individual-level GRM heritability, GCTA GREML/GREML-LDMS and LDAK REML/HE/PCGC
    //
    // Relatedness-matrix construction stays above on the spine so each scientifically distinct matrix is built
    // once and fanned out to every consumer across every domain — this controller and the bivariate
    // relationship controller below both consume matrices built by PREPARE_RELATEDNESS_MATRICES. The
    // controller owns estimator selection, the adaptation of the prepared matrix and headerless phenotype
    // streams into each family's native call shape, and the adjustment-covariate routing only LDAK needs.
    // GCTA and LDAK keep separate native result contracts and are never merged into one heritability table.
    //
    // The two GCTA matrix streams are narrowed to the unary analysis rows here, mirroring the
    // relationship-scoped narrowing the bivariate route below does, so the controller never sees a
    // relationship matrix. An LDAK kinship matrix is only ever requested by a unary heritability method, so
    // that stream is passed as PREPARE_RELATEDNESS_MATRICES emits it.
    ROUTE_GRM_HERITABILITY(
        PREPARE_RELATEDNESS_MATRICES.out.gcta_dense.filter { meta, _grm_files -> !meta.relationship_id },
        PREPARE_RELATEDNESS_MATRICES.out.gcta_ldms.filter { meta, _grm_files, _grm_prefixes -> !meta.relationship_id },
        PREPARE_RELATEDNESS_MATRICES.out.ldak_kinship,
        ch_gcta_phenotypes,
        ch_prepared_adjustment_covariates,
    )

    //
    // SUBWORKFLOW: Pipeline route for MPH REML and REML-LDMS heritability on native MPH matrices
    //
    // A second matrix-backed heritability family. Its matrices are built above on the spine beside the GCTA
    // and LDAK ones and are never interchangeable with them: MPH's layout is its own and both tools read a
    // foreign bundle to completion at exit 0. Its stratified family shares one LD-by-MAF component plan with
    // GCTA's, which is why that plan is now built once on the spine rather than inside either matrix builder.
    //
    // The controller receives five streams because MPH's serializer needs more than the shared headerless
    // seam: the matrix records carry their identity in a tuple position, since a unary analysis row never gets
    // a matrix key and the provenance sidecar needs the plan key and the declared component counts too; and
    // the PLINK 1 bundle is passed because proving the matrix and the genotypes are the same view, and
    // resolving each IID's family identifier, both need the cohort FAM. The matrix streams are narrowed to
    // the unary analysis rows here, mirroring the narrowing the GRM heritability route above receives.
    ROUTE_MPH_HERITABILITY(
        PREPARE_RELATEDNESS_MATRICES.out.mph_dense.filter { meta, _matrix_identity, _grm_files -> !meta.relationship_id },
        PREPARE_RELATEDNESS_MATRICES.out.mph_ldms.filter { meta, _matrix_identity, _grm_files, _grm_prefixes -> !meta.relationship_id },
        PREPARE_COHORT_GENOTYPES.out.plink1_genotypes.filter { meta, _bed, _bim, _fam -> !meta.relationship_id },
        ch_gcta_phenotypes,
        ch_prepared_covariate_tables,
    )

    //
    // SUBWORKFLOW: Pipeline route for LDAK direct-genotype heritability, fast HE and fast PCGC
    //
    // These estimators consume the cohort's PLINK 1 derivative and build no relatedness matrix at all, so they
    // bypass PREPARE_RELATEDNESS_MATRICES entirely: the registry declares no matrix kind for them, and the
    // spine passes the genotype, phenotype and weights streams unmodified. The PLINK 1 stream is narrowed to
    // the unary analysis rows here, mirroring the matrix-stream narrowing above, so the controller never sees
    // a relationship-scoped bundle. The optional LDAK weights Path is narrowed out of the validated row here
    // for the same reason the predictor list is for the association route: reading that row is spine
    // knowledge, and what the weights mean to the estimator is the controller's.
    ROUTE_LDAK_DIRECT_HERITABILITY(
        PREPARE_COHORT_GENOTYPES.out.plink1_genotypes.filter { meta, _bed, _bim, _fam -> !meta.relationship_id },
        ch_gcta_phenotypes,
        ch_analyses.map { meta, _genotype_files, _phenotype, _quant_covariates, _cat_covariates, _kvik_extract, ldak_weights ->
            [meta, ldak_weights ?: []]
        },
    )

    //
    // SUBWORKFLOW: Pipeline route for GCTA bivariate REML and REML-LDMS relationship requests
    //
    // Individual-level relationships are their own domain, disjoint from the summary-statistics pair requests
    // routed below. Dense and LDMS share one controller because they share the relationship definition, the
    // endpoint resolution against the prepared phenotypes and the bivariate trait table built from them. The
    // spine keeps matrix and phenotype preparation; the controller owns relationship
    // de-duplication, declared orientation, preparation reuse and native identity. The matrix
    // streams are narrowed to the relationship-scoped rows here, mirroring the unary narrowing above, so the
    // controller never sees a unary analysis matrix.
    ROUTE_GCTA_BIVARIATE_RELATIONSHIPS(
        ch_relationships,
        PREPARE_RELATIONSHIP_TRAITS.out.pairs,
        PREPARE_RELATEDNESS_MATRICES.out.gcta_dense.filter { meta, _grm_files -> meta.relationship_id },
        PREPARE_RELATEDNESS_MATRICES.out.gcta_ldms.filter { meta, _grm_files, _grm_prefixes -> meta.relationship_id },
    )

    //
    // SUBWORKFLOW: Pipeline route for MPH bivariate REML and REML-LDMS relationship requests
    //
    // The same relationship domain as the GCTA pairs above, over the same prepared pair, but on native MPH
    // matrices, which are a different byte format and are never interchanged with GCTA's. The matrix streams
    // are narrowed to the relationship-scoped rows here, mirroring the unary narrowing the MPH heritability
    // route receives, and the headered covariate stream is the second serialisation of that one prepared pair:
    // MPH names its covariates on the command line, so its serializer reads the header the GCTA estimators
    // reject.
    ROUTE_MPH_BIVARIATE_RELATIONSHIPS(
        ch_relationships,
        PREPARE_RELATIONSHIP_TRAITS.out.pairs,
        PREPARE_RELATIONSHIP_TRAITS.out.named_covariates,
        PREPARE_COHORT_GENOTYPES.out.plink1_genotypes.filter { meta, _bed, _bim, _fam -> meta.relationship_id },
        PREPARE_RELATEDNESS_MATRICES.out.mph_dense.filter { meta, _matrix_identity, _grm_files -> meta.relationship_id },
        PREPARE_RELATEDNESS_MATRICES.out.mph_ldms.filter { meta, _matrix_identity, _grm_files, _grm_prefixes -> meta.relationship_id },
    )

    //
    // SUBWORKFLOW: Pipeline route for GWASLab-standard summary statistics from every origin
    //
    // The single convergence point of the summary-statistics half of the pipeline: it takes the raw
    // association results produced above and the externally supplied sources from the validated manifest, and
    // emits one GWASLab result per summary_statistics_id. The controller owns the internal producer metadata
    // and producer-specific GWASLab mappings. The spine keeps the seam between the association
    // controller above and the fan-out below, and resolves the build-keyed GWASLab resources here because
    // they are pipeline parameters rather than request-owned references.
    def gwaslab_references = getGwaslabReferences()

    ROUTE_CANONICAL_SUMMARY_STATISTICS(
        ROUTE_ASSOCIATION_ANALYSES.out.association_results,
        ch_external_summary_statistics,
        gwaslab_references,
    )

    def ch_selected_meta_requests = ch_meta_requests.filter { meta, _source_ids -> meta.method == 'common_variant_meta_analysis' }
    ROUTE_META_ANALYSIS(ch_selected_meta_requests, ROUTE_CANONICAL_SUMMARY_STATISTICS.out.summary_statistics)
    def ch_summary_statistics = ROUTE_CANONICAL_SUMMARY_STATISTICS.out.summary_statistics
        .mix(ROUTE_META_ANALYSIS.out.summary_statistics)

    //
    // SUBWORKFLOW: Pipeline route for LDAK SumHer and SumCors from GWASLab-standard summary statistics
    //
    // First sibling on the summary-statistics fan-out. Both summary-scale LDAK methods share one
    // controller because they share the GWASLab-to-LDAK preparation and the endpoint resolution that feeds
    // it. The spine selects the route; the controller owns preparation reuse, ordered pair resolution,
    // native-argument and runtime policy, and native outputs. It receives the full validated request tuple so
    // the reference-bundle convention stays request-owned rather than becoming spine knowledge.
    def ch_sumher_requests = ch_unary_requests.filter { meta, _hapmap3_snplist, _reference_ld_scores, _reference_prefix, _regression_weights, _weights_prefix, _tagging_file -> meta.method == 'ldak_sumher' }
    def ch_sumcors_requests = ch_pair_requests.filter { meta, _hapmap3_snplist, _reference_ld_scores, _reference_prefix, _regression_weights, _weights_prefix, _tagging_file -> meta.method == 'ldak_sumcors' }

    ROUTE_LDAK_SUMMARY_ANALYSES(
        ch_sumher_requests,
        ch_sumcors_requests,
        ch_summary_statistics,
    )

    //
    // SUBWORKFLOW: Pipeline route for standalone CBIIT Python 3 LDSC munging, H2 and RG
    //
    // Second sibling on the summary-statistics fan-out. Both summary-scale LDSC methods share one
    // controller because they share the content-addressed munging that feeds them: a GWASLab summary
    // consumed by a unary H2 request and by either side of any number of RG requests is munged exactly once.
    // The spine selects the route; the controller owns the munging reuse identity, endpoint resolution in
    // declared pair order, observed- and liability-scale selection, and native logs. It
    // receives the full validated request tuple so the reference-bundle convention stays request-owned rather
    // than becoming spine knowledge.
    def ch_ldsc_h2_requests = ch_unary_requests.filter { meta, _hapmap3_snplist, _reference_ld_scores, _reference_prefix, _regression_weights, _weights_prefix, _tagging_file -> meta.method == 'ldsc_h2' }
    def ch_ldsc_rg_requests = ch_pair_requests.filter { meta, _hapmap3_snplist, _reference_ld_scores, _reference_prefix, _regression_weights, _weights_prefix, _tagging_file -> meta.method == 'ldsc_rg' }

    ROUTE_LDSC_SUMMARY_ANALYSES(
        ch_ldsc_h2_requests,
        ch_ldsc_rg_requests,
        ch_summary_statistics,
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
    summary_statistics  = ch_summary_statistics // channel: [ val(meta), path(gwaslab_summary_statistics) ]
    multiqc_report      = ROUTE_GWAS_REPORTING.out.report.toList() // channel: [ [ path(report) ] ]
    genotype_views      = ch_genotype_view_records // channel: [ val(cohort_id), path(genotype_view_record) ], one per cohort
    gcta_ldms_artifacts = PREPARE_RELATEDNESS_MATRICES.out.gcta_ldms_artifacts // channel: [ val(matrix_meta), path(grm_files), val(grm_prefixes) ], one per base key
    mph_dense_artifacts = PREPARE_RELATEDNESS_MATRICES.out.mph_dense_artifacts // channel: [ val(matrix_meta), path(grm_files) ], one per base key
    mph_ldms_artifacts  = PREPARE_RELATEDNESS_MATRICES.out.mph_ldms_artifacts // channel: [ val(matrix_meta), path(grm_files), val(grm_prefixes) ], one per base key
    ldms_plan_artifacts = PREPARE_RELATEDNESS_MATRICES.out.ldms_plan_artifacts // channel: [ val(plan_meta), path(ld_scores), path(strata_manifest), [ path(snp_group_file), ... ] ], one per plan key
}
