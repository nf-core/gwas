/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// MODULE: Local to the pipeline
include { CANONICALISE_SUMMARY_STATISTICS                     } from '../modules/local/canonicalise_summary_statistics/main'
include { GWASLAB_HARMONIZE                                   } from '../modules/local/gwaslab/harmonize/main'
include { GCTA_FASTGWA                                        } from '../modules/local/gcta/fastgwa/main'
include { LDSC_H2 as LDSC_H2_LIABILITY                        } from '../modules/local/ldsc/h2/main'
include { LDSC_H2 as LDSC_H2_OBSERVED                         } from '../modules/local/ldsc/h2/main'
include { LDSC_MUNGESUMSTATS                                  } from '../modules/local/ldsc/mungesumstats/main'
include { LDSC_RG as LDSC_RG_LIABILITY                        } from '../modules/local/ldsc/rg/main'
include { LDSC_RG as LDSC_RG_OBSERVED                         } from '../modules/local/ldsc/rg/main'
include { LDAK_SUMCORS                                        } from '../modules/local/ldak/sumcors/main'
include { LDAK_SUMHER                                         } from '../modules/local/ldak/sumher/main'
include { NORMALISE_LDAK_SUMCORS                              } from '../modules/local/normalise_ldak_sumcors/main'
include { NORMALISE_LDAK_SUMHER                               } from '../modules/local/normalise_ldak_sumher/main'
include { NORMALISE_PHENOTYPES                                } from '../modules/local/normalise_phenotypes/main'
include { NORMALISE_GCTA_BIVARIATE                            } from '../modules/local/normalise_gcta_bivariate/main'
include { NORMALISE_LDSC                                      } from '../modules/local/normalise_ldsc/main'
include { PLINK2_GLM                                          } from '../modules/local/plink2/glm/main'
include { PREPARE_BIVARIATE_TRAITS                            } from '../modules/local/prepare_bivariate_traits/main'
include { PREPARE_LDAK_SUMMARY_STATISTICS                     } from '../modules/local/prepare_ldak_summary_statistics/main'

// MODULE: Installed directly from nf-core/modules
include { GCTA_BIVARIATEREML                                  } from '../modules/nf-core/gcta/bivariatereml/main'
include { GCTA_BIVARIATEREMLLDMS                              } from '../modules/nf-core/gcta/bivariateremlldms/main'
include { MULTIQC                                             } from '../modules/nf-core/multiqc/main'

// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
include { GRM_HERITABILITY_GCTA                               } from '../subworkflows/local/grm_heritability_gcta'
include { GRM_HERITABILITY_LDAK as GRM_HERITABILITY_LDAK_HE   } from '../subworkflows/local/grm_heritability_ldak'
include { GRM_HERITABILITY_LDAK as GRM_HERITABILITY_LDAK_PCGC } from '../subworkflows/local/grm_heritability_ldak'
include { GRM_HERITABILITY_LDAK as GRM_HERITABILITY_LDAK_REML } from '../subworkflows/local/grm_heritability_ldak'
include { PREPARE_COHORT_GENOTYPES                            } from '../subworkflows/local/prepare_cohort_genotypes'
include { PREPARE_RELATEDNESS_MATRICES                        } from '../subworkflows/local/prepare_relatedness_matrices'
include { ROUTE_LDAK_KVIK_ASSOCIATIONS                        } from '../subworkflows/local/route_ldak_kvik_associations'
include { ROUTE_REGENIE_ASSOCIATIONS                          } from '../subworkflows/local/route_regenie_associations'
include { getAssociationColumnMappingJson                     } from '../subworkflows/local/validate_gwas_input'
include { getInternalSummaryMetadata                          } from '../subworkflows/local/validate_gwas_input'
include { getGwaslabReferences                                } from '../subworkflows/local/utils_nfcore_gwas_pipeline'
include { analysisPlanJson                                    } from '../subworkflows/local/utils_nfcore_gwas_pipeline'
include { digestFileBytes                                     } from '../subworkflows/local/utils_nfcore_gwas_pipeline'
include { digestIdentityText                                  } from '../subworkflows/local/utils_nfcore_gwas_pipeline'
include { methodsDescriptionText                              } from '../subworkflows/local/utils_nfcore_gwas_pipeline'

// SUBWORKFLOW: Consisting entirely of nf-core/modules
include { paramsSummaryMultiqc                                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'

// PLUGIN
include { paramsSummaryMap                                    } from 'plugin/nf-schema'

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
    ch_unary_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ]
    ch_pair_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ]
    multiqc_config // channel: val(multiqc_config)
    multiqc_logo // channel: val(multiqc_logo)
    multiqc_methods_description // channel: val(multiqc_methods_description)
    outdir // channel: val(outdir)

    main:

    def ch_multiqc_files = channel.empty()
    def ch_analysis_metadata = ch_analyses
        .map { meta, _genotype_files, _phenotype, _quant_covariates, _cat_covariates, _kvik_extract, _ldak_weights -> meta }
        .collect()
    def ch_method_metadata = ch_analyses
        .map { meta, _genotype_files, _phenotype, _quant_covariates, _cat_covariates, _kvik_extract, _ldak_weights -> [domain: 'analysis', meta: meta] }
        .mix(
            ch_relationships.map { meta, _genotype_files, _pair_quant_covariates, _pair_cat_covariates -> [domain: 'pairwise', meta: meta] }
        )
        .mix(
            ch_unary_requests.map { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file -> [domain: 'summary_unary', meta: meta] }
        )
        .mix(
            ch_pair_requests.map { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file -> [domain: 'pairwise', meta: meta] }
        )
        .collect()

    // One element per analysis unit carrying the genotype files it declared. A pairwise LDMS request is
    // also a consumer of the cohort's lazy PLINK 1 derivative, so it enters this request stream without
    // inheriting either endpoint's unary method settings. Cohort preparation collapses every request to
    // the distinct cohort before conversion.
    def ch_genotype_requests = ch_analyses.map { meta, genotype_files, _phenotype, _quant_covariates, _cat_covariates, _kvik_extract, _ldak_weights ->
        [meta, genotype_files]
    }
    ch_genotype_requests = ch_genotype_requests.mix(
        ch_relationships
            .filter { meta, _genotype_files, _pair_quant_covariates, _pair_cat_covariates -> meta.matrix_kind == 'gcta_ldms' }
            .map { meta, genotype_files, _pair_quant_covariates, _pair_cat_covariates -> [meta, genotype_files] }
    )

    // Matrix preparation additionally receives the optional LDAK weights Path. It derives identity from the
    // bytes before request deduplication and keeps the Path outside matrix metadata and the published key.
    def ch_relatedness_analyses = ch_analyses.map { meta, genotype_files, _phenotype, _quant_covariates, _cat_covariates, _kvik_extract, ldak_weights ->
        [meta, genotype_files, ldak_weights ?: []]
    }
    ch_relatedness_analyses = ch_relatedness_analyses.mix(
        ch_relationships.map { meta, genotype_files, _pair_quant_covariates, _pair_cat_covariates -> [meta, genotype_files, []] }
    )

    //
    // SUBWORKFLOW: Prepare each distinct cohort's genotypes once into the canonical PLINK 2 bundle
    //
    PREPARE_COHORT_GENOTYPES(ch_genotype_requests)

    //
    // SUBWORKFLOW: Build each distinct relatedness matrix once and fan it out per analysis unit
    //
    PREPARE_RELATEDNESS_MATRICES(
        ch_relatedness_analyses,
        PREPARE_COHORT_GENOTYPES.out.cohort_genotypes,
        PREPARE_COHORT_GENOTYPES.out.plink1_genotypes,
        params.gcta_grm_parts,
    )

    //
    // MODULE: Normalise each analysis unit's phenotype and covariates into the canonical layout
    //
    NORMALISE_PHENOTYPES(
        ch_analyses.map { meta, _genotype_files, phenotype, quant_covariates, cat_covariates, _kvik_extract, _ldak_weights ->
            [meta, phenotype, quant_covariates, cat_covariates]
        }
    )

    // The genotype bundle, the normalised phenotype and the merged covariate file, one element per
    // analysis unit. `join` is correct here where `combine` was correct at the cohort seam: all three
    // channels are keyed one-to-one on the analysis meta, so a missing or duplicated key is a defect
    // and the strict form is what says so. The covariate file is optional, so it joins with
    // `remainder: true` and arrives as `null` for a row that supplied none.
    def ch_analysis_inputs = PREPARE_COHORT_GENOTYPES.out.genotypes
        .filter { meta, _pgen, _psam, _pvar -> !meta.relationship_id }
        .join(NORMALISE_PHENOTYPES.out.phenotype, failOnMismatch: true, failOnDuplicate: true)
        .join(NORMALISE_PHENOTYPES.out.covariates, remainder: true)

    // GCTA and LDAK reject a header row. fastGWA, GREML and LDAK REML therefore consume the headerless
    // phenotype and covariate serialisations. Optional covariates are represented by [], which stages nothing.
    def ch_gcta_phenotypes = NORMALISE_PHENOTYPES.out.phenotype_headerless
        .join(NORMALISE_PHENOTYPES.out.quant_covariates_headerless, remainder: true)
        .join(NORMALISE_PHENOTYPES.out.cat_covariates_headerless, remainder: true)
        .map { meta, phenotype, quant_covariates, cat_covariates ->
            [meta, phenotype, quant_covariates ?: [], cat_covariates ?: []]
        }

    // Relationships own their orientation and covariates. Collapse the per-method request fan-out to one
    // relationship definition, resolve each endpoint against the canonical unary phenotype stream, and
    // construct one ordered full-union two-trait table. `combine` is deliberate at the endpoint seams: one
    // analysis may be reused by several relationships. The prepared artifact is fanned back out by
    // relationship ID only after construction, so selecting dense and LDMS does not duplicate it.
    def ch_relationship_definitions = ch_relationships
        .map { meta, genotype_files, pair_quant_covariates, pair_cat_covariates ->
            def relationship_meta = meta + [
                id: meta.relationship_id,
                request_id: meta.relationship_id,
            ]
            [meta.relationship_id, relationship_meta, genotype_files, pair_quant_covariates, pair_cat_covariates]
        }
        .unique { relationship_id, _meta, _genotype_files, _pair_quant_covariates, _pair_cat_covariates -> relationship_id }

    def ch_left_pair_phenotypes = ch_relationship_definitions
        .map { relationship_id, meta, _genotype_files, pair_quant_covariates, pair_cat_covariates ->
            [meta.left_analysis_id, relationship_id, meta, pair_quant_covariates ?: [], pair_cat_covariates ?: []]
        }
        .combine(
            NORMALISE_PHENOTYPES.out.phenotype_headerless.map { meta, phenotype -> [meta.id, phenotype] },
            by: 0
        )
        .map { _analysis_id, relationship_id, meta, pair_quant_covariates, pair_cat_covariates, phenotype ->
            [relationship_id, meta, phenotype, pair_quant_covariates, pair_cat_covariates]
        }

    def ch_right_pair_phenotypes = ch_relationship_definitions
        .map { relationship_id, meta, _genotype_files, _pair_quant_covariates, _pair_cat_covariates -> [meta.right_analysis_id, relationship_id] }
        .combine(
            NORMALISE_PHENOTYPES.out.phenotype_headerless.map { meta, phenotype -> [meta.id, phenotype] },
            by: 0
        )
        .map { _analysis_id, relationship_id, phenotype -> [relationship_id, phenotype] }

    def ch_pair_trait_inputs = ch_left_pair_phenotypes
        .join(ch_right_pair_phenotypes, failOnDuplicate: true, failOnMismatch: true)
        .map { _relationship_id, meta, left_phenotype, pair_quant_covariates, pair_cat_covariates, right_phenotype ->
            [meta, left_phenotype, right_phenotype, pair_quant_covariates, pair_cat_covariates]
        }

    PREPARE_BIVARIATE_TRAITS(ch_pair_trait_inputs)

    //
    // MODULE: PLINK 2 --glm association
    //
    // `multiMap` rather than three `map`s of the same channel, so the three inputs cannot drift out
    // of lockstep. A row that supplied no covariates passes `[]`, which stages nothing: the module's
    // covariate argument is a ternary on a `path` inside a tuple, and no placeholder file is written.
    def ch_glm_input = ch_analysis_inputs
        .filter { meta, _pgen, _psam, _pvar, _phenotype, _covariates -> 'plink2' in meta.association_methods }
        .multiMap { meta, pgen, psam, pvar, phenotype, covariates ->
            genotypes: [meta, pgen, psam, pvar]
            phenotype: [meta, phenotype]
            covariates: [meta, covariates ?: []]
        }

    PLINK2_GLM(
        ch_glm_input.genotypes,
        ch_glm_input.phenotype,
        ch_glm_input.covariates,
    )

    //
    // PIPELINE ROUTE: REGENIE association with shared Step 1 predictions
    //
    // The local route owns nf-core/gwas scientific identity, cross-analysis fit reuse and output
    // attribution. Upstream-ready REGENIE components remain unaware of the relational input contract.
    def ch_regenie_analyses = ch_analysis_inputs.filter { meta, _pgen, _psam, _pvar, _phenotype, _covariates -> 'regenie' in meta.association_methods }

    ROUTE_REGENIE_ASSOCIATIONS(
        ch_regenie_analyses,
        params.regenie_step2_bsize,
        params.regenie_step1_mode,
        params.regenie_step1_jobs,
    )

    //
    // PIPELINE ROUTE: LDAK-KVIK association with shared Step 1 predictions
    //
    // LDAK consumes the headerless phenotype serialisation and keeps quantitative and categorical
    // covariates separate. Fold both optional covariate streams onto the total phenotype stream so
    // that an absent file is represented by `[]` and stages nothing.
    def ch_kvik_phenotypes = NORMALISE_PHENOTYPES.out.phenotype_headerless
        .join(NORMALISE_PHENOTYPES.out.quant_covariates_headerless, remainder: true)
        .join(NORMALISE_PHENOTYPES.out.cat_covariates_headerless, remainder: true)
        .filter { meta, _phenotype, _quant_covariates, _cat_covariates -> 'ldak_kvik' in meta.association_methods }
        .map { meta, phenotype, quant_covariates, cat_covariates ->
            [meta.id, meta, phenotype, quant_covariates ?: [], cat_covariates ?: []]
        }

    // The stageable predictor resource and its validated policy come from the relational LDAK family map.
    // Both remain explicit tuple members so Nextflow stages the file, while their content identity is folded
    // into the Step 1 reuse key below.
    def ch_kvik_extract_policy = ch_analyses
        .filter { meta, _genotype_files, _phenotype, _quant_covariates, _cat_covariates, _kvik_extract, _ldak_weights -> 'ldak_kvik' in meta.association_methods }
        .map { meta, _genotype_files, _phenotype, _quant_covariates, _cat_covariates, kvik_extract, _ldak_weights ->
            [meta.id, meta, kvik_extract ?: [], meta.method_options.ldak.kvik_step1_subset]
        }

    def ch_kvik_genotypes = PREPARE_COHORT_GENOTYPES.out.plink1_genotypes.filter { meta, _bed, _bim, _fam -> 'ldak_kvik' in (meta.association_methods ?: []) }

    ROUTE_LDAK_KVIK_ASSOCIATIONS(
        ch_kvik_genotypes,
        ch_kvik_phenotypes.map { _analysis_id, meta, phenotype, quant_covariates, cat_covariates -> [meta, phenotype, quant_covariates, cat_covariates] },
        ch_kvik_extract_policy.map { _analysis_id, meta, kvik_extract, subset_policy -> [meta, kvik_extract, subset_policy] },
    )

    //
    // MODULE: GCTA fastGWA-MLM association
    //
    // This is deliberately inline: a composition wrapping one module is not an nf-core subworkflow.
    // The module chooses --fastGWA-mlm or --fastGWA-mlm-binary from the boolean phenotype input;
    // conf/modules/gcta.config supplies no arbitrary ext.args, so plain --fastGWA-lr is unreachable.
    def ch_fastgwa_genotypes = PREPARE_COHORT_GENOTYPES.out.genotypes.filter { meta, _pgen, _psam, _pvar -> 'gcta_fastgwa' in meta.association_methods }
    def ch_fastgwa_phenotypes = ch_gcta_phenotypes.filter { meta, _phenotype, _quant_covariates, _cat_covariates -> 'gcta_fastgwa' in meta.association_methods }

    def ch_fastgwa_invocations = ch_fastgwa_genotypes
        .map { meta, pgen, psam, pvar -> [meta.id, [meta, pgen, pvar, psam]] }
        .join(
            ch_fastgwa_phenotypes.map { meta, phenotype, _quant_covariates, _cat_covariates -> [meta.id, [meta, phenotype, meta.is_binary]] },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_fastgwa_phenotypes.map { meta, _phenotype, quant_covariates, _cat_covariates -> [meta.id, [meta, quant_covariates]] },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_fastgwa_phenotypes.map { meta, _phenotype, _quant_covariates, cat_covariates -> [meta.id, [meta, cat_covariates]] },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            PREPARE_RELATEDNESS_MATRICES.out.gcta_sparse.map { meta, sparse_grm_files -> [meta.id, [meta, sparse_grm_files]] },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .multiMap { _analysis_id, genotypes, pheno, qcovar, covar, sparse_grm ->
            genotypes: genotypes
            pheno: pheno
            qcovar: qcovar
            covar: covar
            sparse_grm: sparse_grm
        }

    GCTA_FASTGWA(
        ch_fastgwa_invocations.genotypes,
        ch_fastgwa_invocations.pheno,
        ch_fastgwa_invocations.qcovar,
        ch_fastgwa_invocations.covar,
        ch_fastgwa_invocations.sparse_grm,
    )

    //
    // MODULE: GWASLab harmonisation of every association result
    //
    // One record per analysis per association method actually exercised. Each route contributes an
    // adapter that names its method on the meta map and normalises whatever emissions the programme
    // splits its results across; everything downstream is method-agnostic. `meta.id` stays the analysis
    // identifier — the method is a separate key, because the analysis is what the published summary
    // statistics directory is keyed by and the method is what distinguishes the files inside it.
    //
    // PLINK 2 splits its result across four optional emissions, one per regression it may have fitted,
    // and exactly one of them is populated for a given analysis, so the four are mixed back into one.
    def ch_association_results = channel.empty()

    def ch_plink2_results = PLINK2_GLM.out.linear
    ch_plink2_results = ch_plink2_results.mix(PLINK2_GLM.out.logistic)
    ch_plink2_results = ch_plink2_results.mix(PLINK2_GLM.out.logistic_hybrid)
    ch_plink2_results = ch_plink2_results.mix(PLINK2_GLM.out.firth)

    ch_association_results = ch_association_results.mix(
        ch_plink2_results.map { meta, sumstats -> [meta + [method: 'plink2'], sumstats] }
    )
    ch_association_results = ch_association_results.mix(
        ROUTE_REGENIE_ASSOCIATIONS.out.results.map { meta, sumstats -> [meta + [method: 'regenie'], sumstats] }
    )
    ch_association_results = ch_association_results.mix(
        ROUTE_LDAK_KVIK_ASSOCIATIONS.out.harmonisation_input.map { meta, sumstats -> [meta + [method: 'ldak_kvik'], sumstats] }
    )
    ch_association_results = ch_association_results.mix(
        GCTA_FASTGWA.out.results.map { meta, sumstats -> [meta + [method: 'gcta_fastgwa'], sumstats] }
    )

    // Internal association results and external raw inputs converge before GWASLab. The producer-specific
    // internal mappings remain explicit and unchanged; an external row supplies a named GWASLab format.
    // Already-canonical external inputs bypass GWASLab and enter only the canonical contract validator.
    def ch_harmonise_records = ch_association_results
        .map { meta, source ->
            def summary_meta = getInternalSummaryMetadata(meta, meta.method) + [
                method: meta.method,
                source_name: source.name,
                gwaslab_input_format: getAssociationColumnMappingJson(meta.method, meta.is_binary),
            ]
            [summary_meta, source]
        }
        .mix(
            ch_external_summary_statistics.filter { meta, _source -> meta.source_mode == 'raw' }.map { meta, source -> [meta + [method: 'external', gwaslab_input_format: meta.source_format], source] }
        )

    // The optional GWASLab resources remain build-keyed pipeline parameters. This is independent from the
    // request-owned LDSC/LDAK reference catalog and does not infer a scientific analysis reference.
    def gwaslab_references = getGwaslabReferences()
    def ch_harmonise_input = ch_harmonise_records.multiMap { meta, source ->
        def references = gwaslab_references[meta.build]
        sumstats: [meta, source, meta.gwaslab_input_format, meta.build]
        reference_fasta: [[id: meta.build], references.fasta, references.fasta_index]
        rsid_reference: [[id: meta.build], references.rsid_vcf, references.rsid_vcf_index]
        strand_reference: [[id: meta.build], references.strand_vcf, references.strand_vcf_index]
    }

    GWASLAB_HARMONIZE(
        ch_harmonise_input.sumstats,
        ch_harmonise_input.reference_fasta,
        ch_harmonise_input.rsid_reference,
        ch_harmonise_input.strand_reference,
    )

    def ch_harmonise_sources = ch_harmonise_records.map { meta, source -> [meta.summary_statistics_id, source] }
    def ch_canonical_candidates = GWASLAB_HARMONIZE.out.sumstats
        .map { meta, candidate -> [meta.summary_statistics_id, meta, candidate] }
        .join(ch_harmonise_sources, failOnDuplicate: true, failOnMismatch: true)
        .map { _summary_statistics_id, meta, candidate, source -> [meta, candidate, source] }
        .mix(
            ch_external_summary_statistics.filter { meta, _source -> meta.source_mode == 'canonical' }.map { meta, source -> [meta, source, source] }
        )

    CANONICALISE_SUMMARY_STATISTICS(ch_canonical_candidates)

    //
    // PIPELINE ROUTES: LDAK SumHer and SumCors from canonical summary statistics
    //
    // Adapt each distinct summary once, regardless of how many unary, pairwise or named sensitivity
    // requests consume it. The adapter owns only the deterministic canonical-to-LDAK column transform;
    // each request retains its own tagging reference, effective native arguments and publication identity.
    def ch_ldak_requested_summary_ids = ch_unary_requests
        .filter { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file -> meta.method == 'ldak_sumher' }
        .map { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file -> [meta.summary_statistics_id] }
        .mix(
            ch_pair_requests.filter { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file -> meta.method == 'ldak_sumcors' }.flatMap { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file ->
                [[meta.left_summary_statistics_id], [meta.right_summary_statistics_id]]
            }
        )
        .unique()

    def ch_ldak_canonical_summaries = ch_ldak_requested_summary_ids
        .combine(
            CANONICALISE_SUMMARY_STATISTICS.out.summary_statistics.map { meta, summary_statistics -> [meta.summary_statistics_id, meta, summary_statistics] },
            by: 0
        )
        .map { _summary_statistics_id, meta, summary_statistics -> [meta, summary_statistics] }

    PREPARE_LDAK_SUMMARY_STATISTICS(ch_ldak_canonical_summaries)

    def ch_prepared_ldak_summaries = PREPARE_LDAK_SUMMARY_STATISTICS.out.summary_statistics
        .map { meta, summary_statistics -> [meta.summary_statistics_id, meta, summary_statistics] }
        .join(
            PREPARE_LDAK_SUMMARY_STATISTICS.out.preparation.map { meta, preparation -> [meta.summary_statistics_id, preparation] },
            failOnDuplicate: true,
            failOnMismatch: true,
        )

    def ch_sumher_invocations = ch_unary_requests
        .filter { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file -> meta.method == 'ldak_sumher' }
        .map { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, tagging_file -> [meta.summary_statistics_id, meta, tagging_file] }
        .combine(ch_prepared_ldak_summaries, by: 0)
        .multiMap { summary_statistics_id, request_meta, tagging_file, _summary_meta, summary_statistics, preparation ->
            def route_meta = request_meta + [
                id: summary_statistics_id,
                effective_native_args: getLdakSummaryArguments(request_meta),
                native_runtime: getLdakSummaryRuntime(),
            ]
            summary: [route_meta, summary_statistics]
            tagging: [[id: request_meta.reference_bundle_id], tagging_file]
            preparation: [request_meta.request_id, preparation]
        }

    LDAK_SUMHER(
        ch_sumher_invocations.summary,
        ch_sumher_invocations.tagging,
    )

    def ch_sumher_native_results = LDAK_SUMHER.out.hers
        .map { meta, hers -> [meta.request_id, meta, hers] }
        .join(LDAK_SUMHER.out.extra.map { meta, extra -> [meta.request_id, extra] }, failOnDuplicate: true, failOnMismatch: true)
        .join(LDAK_SUMHER.out.overlap.map { meta, overlap -> [meta.request_id, overlap] }, failOnDuplicate: true, failOnMismatch: true)
        .join(LDAK_SUMHER.out.log.map { meta, ldak_log -> [meta.request_id, ldak_log] }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_sumher_invocations.preparation, failOnDuplicate: true, failOnMismatch: true)
        .join(
            LDAK_SUMHER.out.hers_liability.map { meta, hers_liability -> [meta.request_id, hers_liability] },
            remainder: true,
            failOnDuplicate: true,
        )
        .map { _request_id, meta, hers, extra, overlap, ldak_log, preparation, hers_liability ->
            [meta, hers, extra, overlap, ldak_log, preparation, hers_liability ?: []]
        }

    NORMALISE_LDAK_SUMHER(ch_sumher_native_results)

    def ch_sumcors_left = ch_pair_requests
        .filter { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file -> meta.method == 'ldak_sumcors' }
        .map { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, tagging_file -> [meta.left_summary_statistics_id, meta, tagging_file] }
        .combine(ch_prepared_ldak_summaries, by: 0)
        .map { _left_summary_statistics_id, request_meta, tagging_file, _left_meta, left_summary_statistics, left_preparation ->
            [request_meta.right_summary_statistics_id, request_meta, tagging_file, left_summary_statistics, left_preparation]
        }

    def ch_sumcors_invocations = ch_sumcors_left
        .combine(ch_prepared_ldak_summaries, by: 0)
        .multiMap { _right_summary_statistics_id, request_meta, tagging_file, left_summary_statistics, left_preparation, right_meta, right_summary_statistics, right_preparation ->
            def route_meta = request_meta + [
                id: request_meta.left_summary_statistics_id,
                effective_native_args: getLdakSummaryArguments(request_meta),
                native_runtime: getLdakSummaryRuntime(),
            ]
            left: [route_meta, left_summary_statistics]
            right: [right_meta, right_summary_statistics]
            tagging: [[id: request_meta.reference_bundle_id], tagging_file]
            preparation: [request_meta.request_id, left_preparation, right_preparation]
        }

    LDAK_SUMCORS(
        ch_sumcors_invocations.left,
        ch_sumcors_invocations.right,
        ch_sumcors_invocations.tagging,
    )

    def ch_sumcors_native_results = LDAK_SUMCORS.out.correlations
        .map { meta, _meta2, correlations -> [meta.request_id, meta, correlations] }
        .join(LDAK_SUMCORS.out.correlations_full.map { meta, _meta2, correlations_full -> [meta.request_id, correlations_full] }, failOnDuplicate: true, failOnMismatch: true)
        .join(LDAK_SUMCORS.out.overlap.map { meta, _meta2, overlap -> [meta.request_id, overlap] }, failOnDuplicate: true, failOnMismatch: true)
        .join(LDAK_SUMCORS.out.log.map { meta, _meta2, ldak_log -> [meta.request_id, ldak_log] }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_sumcors_invocations.preparation, failOnDuplicate: true, failOnMismatch: true)
        .join(
            LDAK_SUMCORS.out.correlations_liability.map { meta, _meta2, correlations_liability -> [meta.request_id, correlations_liability] },
            remainder: true,
            failOnDuplicate: true,
        )
        .map { _request_id, meta, correlations, correlations_full, overlap, ldak_log, left_preparation, right_preparation, correlations_liability ->
            [meta, correlations, correlations_full, overlap, ldak_log, left_preparation, right_preparation, correlations_liability ?: []]
        }

    NORMALISE_LDAK_SUMCORS(ch_sumcors_native_results)

    // PIPELINE ROUTE: standalone CBIIT Python 3 LDSC munging, H2 and RG
    //
    // Munging belongs to a canonical summary plus the exact HapMap3 allele-universe bytes, not to a
    // downstream request or its regression reference/weights. A content-derived key therefore lets unary,
    // pairwise, primary and named sensitivity requests reuse the same expensive preparation without making
    // ancestry, bundle names, H2/RG native arguments or output identity part of that derivation.
    def ldsc_munging_key = { summary_statistics_id, hapmap3_snplist ->
        digestIdentityText([
            'adapter=nfcore_gwas_canonical_v1_to_ldsc_sumstats_v1',
            "summary_statistics_id=${summary_statistics_id}",
            "hapmap3_sha256=${digestFileBytes(hapmap3_snplist)}",
        ].join('\n'))
    }

    def ch_ldsc_munging_requests = ch_unary_requests
        .filter { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file -> meta.method == 'ldsc_h2' }
        .map { meta, hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file ->
            def key = ldsc_munging_key(meta.summary_statistics_id, hapmap3_snplist)
            [meta.summary_statistics_id, key, hapmap3_snplist]
        }
        .mix(
            ch_pair_requests
                .filter { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file -> meta.method == 'ldsc_rg' }
                .flatMap { meta, hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file ->
                    [meta.left_summary_statistics_id, meta.right_summary_statistics_id].collect { summary_statistics_id ->
                        [summary_statistics_id, ldsc_munging_key(summary_statistics_id, hapmap3_snplist), hapmap3_snplist]
                    }
                }
        )
        .unique { _summary_statistics_id, key, _hapmap3_snplist -> key }

    def ch_canonical_by_summary_id = CANONICALISE_SUMMARY_STATISTICS.out.summary_statistics
        .map { meta, canonical_summary_statistics -> [meta.summary_statistics_id, meta, canonical_summary_statistics] }

    def ch_ldsc_munging_invocations = ch_ldsc_munging_requests
        .combine(ch_canonical_by_summary_id, by: 0)
        .multiMap { summary_statistics_id, key, hapmap3_snplist, summary_meta, canonical_summary_statistics ->
            def munging_meta = summary_meta + [
                id: key,
                munging_key: key,
                summary_statistics_id: summary_statistics_id,
                hapmap3_sha256: digestFileBytes(hapmap3_snplist),
                munging_adapter_contract: 'nfcore_gwas_canonical_v1_to_ldsc_sumstats_v1',
            ]
            sumstats: [munging_meta, canonical_summary_statistics]
            merge_alleles: [[id: key], hapmap3_snplist]
        }

    LDSC_MUNGESUMSTATS(
        ch_ldsc_munging_invocations.sumstats,
        ch_ldsc_munging_invocations.merge_alleles,
    )

    def ch_ldsc_munged = LDSC_MUNGESUMSTATS.out.munged_sumstats
        .map { meta, munged_sumstats -> [meta.munging_key, meta, munged_sumstats] }
        .join(
            LDSC_MUNGESUMSTATS.out.log.map { meta, munging_log -> [meta.munging_key, munging_log] },
            failOnDuplicate: true,
            failOnMismatch: true,
        )

    // Unary request identity and request-owned LD/weight resources are joined only after munging. Observed
    // scale is always retained. A second native invocation is made only when a binary endpoint declares both
    // population and sample prevalence, because native LDSC emits liability rather than observed H2 when
    // those values are supplied.
    def ch_ldsc_h2_requests = ch_unary_requests
        .filter { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file -> meta.method == 'ldsc_h2' }
        .map { meta, hapmap3_snplist, reference_ld_scores, regression_weights, _tagging_file ->
            [ldsc_munging_key(meta.summary_statistics_id, hapmap3_snplist), meta, reference_ld_scores, regression_weights]
        }
        .combine(ch_ldsc_munged, by: 0)

    def ch_ldsc_h2_observed = ch_ldsc_h2_requests.multiMap { key, meta, reference_ld_scores, regression_weights, munging_meta, munged_sumstats, munging_log ->
        def route_meta = meta + [munging_keys: [key], native_scale: 'observed']
        sumstats: [route_meta, munged_sumstats]
        reference_ld_scores: [[id: meta.reference_bundle_id], reference_ld_scores]
        regression_weights: [[id: meta.reference_bundle_id], regression_weights]
        munging_log: [meta.request_id, munging_log]
    }

    LDSC_H2_OBSERVED(
        ch_ldsc_h2_observed.sumstats,
        ch_ldsc_h2_observed.reference_ld_scores,
        ch_ldsc_h2_observed.regression_weights,
    )

    def ch_ldsc_h2_liability = ch_ldsc_h2_requests
        .filter { _key, meta, _reference_ld_scores, _regression_weights, _munging_meta, _munged_sumstats, _munging_log ->
            meta.is_binary && meta.population_prevalence != null && meta.sample_prevalence != null
        }
        .multiMap { key, meta, reference_ld_scores, regression_weights, munging_meta, munged_sumstats, munging_log ->
            def route_meta = meta + [
                munging_keys: [key],
                native_scale: 'liability',
                effective_population_prevalence: [meta.population_prevalence],
                effective_sample_prevalence: [meta.sample_prevalence],
            ]
            sumstats: [route_meta, munged_sumstats]
            reference_ld_scores: [[id: meta.reference_bundle_id], reference_ld_scores]
            regression_weights: [[id: meta.reference_bundle_id], regression_weights]
        }

    LDSC_H2_LIABILITY(
        ch_ldsc_h2_liability.sumstats,
        ch_ldsc_h2_liability.reference_ld_scores,
        ch_ldsc_h2_liability.regression_weights,
    )

    // Pair requests preserve declared left/right order. Both endpoint munging keys are resolved against the
    // one HapMap3 resource selected by this request, then the request-owned LD-score and regression-weight
    // directories are passed unchanged to RG.
    def ch_ldsc_rg_left = ch_pair_requests
        .filter { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file -> meta.method == 'ldsc_rg' }
        .map { meta, hapmap3_snplist, reference_ld_scores, regression_weights, _tagging_file ->
            def left_key = ldsc_munging_key(meta.left_summary_statistics_id, hapmap3_snplist)
            [left_key, meta, hapmap3_snplist, reference_ld_scores, regression_weights]
        }
        .combine(ch_ldsc_munged, by: 0)
        .map { left_key, meta, hapmap3_snplist, reference_ld_scores, regression_weights, left_munging_meta, left_sumstats, left_munging_log ->
            def right_key = ldsc_munging_key(meta.right_summary_statistics_id, hapmap3_snplist)
            [right_key, left_key, meta, reference_ld_scores, regression_weights, left_sumstats, left_munging_log]
        }

    def ch_ldsc_rg_requests = ch_ldsc_rg_left
        .combine(ch_ldsc_munged, by: 0)
        .map { right_key, left_key, meta, reference_ld_scores, regression_weights, left_sumstats, left_munging_log, right_munging_meta, right_sumstats, right_munging_log ->
            [meta, left_key, right_key, reference_ld_scores, regression_weights, left_sumstats, right_sumstats, left_munging_log, right_munging_log]
        }

    def ch_ldsc_rg_observed = ch_ldsc_rg_requests.multiMap { meta, left_key, right_key, reference_ld_scores, regression_weights, left_sumstats, right_sumstats, left_munging_log, right_munging_log ->
        def route_meta = meta + [munging_keys: [left_key, right_key], native_scale: 'observed']
        sumstats: [route_meta, left_sumstats, right_sumstats]
        reference_ld_scores: [[id: meta.reference_bundle_id], reference_ld_scores]
        regression_weights: [[id: meta.reference_bundle_id], regression_weights]
        munging_logs: [meta.request_id, [left_munging_log, right_munging_log]]
    }

    LDSC_RG_OBSERVED(
        ch_ldsc_rg_observed.sumstats,
        ch_ldsc_rg_observed.reference_ld_scores,
        ch_ldsc_rg_observed.regression_weights,
    )

    def ch_ldsc_rg_liability = ch_ldsc_rg_requests
        .filter { meta, _left_key, _right_key, _reference_ld_scores, _regression_weights, _left_sumstats, _right_sumstats, _left_munging_log, _right_munging_log ->
            def has_binary = meta.left_is_binary || meta.right_is_binary
            def complete = [
                [binary: meta.left_is_binary, population: meta.left_population_prevalence, sample: meta.left_sample_prevalence],
                [binary: meta.right_is_binary, population: meta.right_population_prevalence, sample: meta.right_sample_prevalence],
            ].every { endpoint -> !endpoint.binary || (endpoint.population != null && endpoint.sample != null) }
            has_binary && complete
        }
        .multiMap { meta, left_key, right_key, reference_ld_scores, regression_weights, left_sumstats, right_sumstats, left_munging_log, right_munging_log ->
            def population = [
                meta.left_is_binary ? meta.left_population_prevalence : 'nan',
                meta.right_is_binary ? meta.right_population_prevalence : 'nan',
            ]
            def sample = [
                meta.left_is_binary ? meta.left_sample_prevalence : 'nan',
                meta.right_is_binary ? meta.right_sample_prevalence : 'nan',
            ]
            def route_meta = meta + [
                munging_keys: [left_key, right_key],
                native_scale: 'liability',
                effective_population_prevalence: population,
                effective_sample_prevalence: sample,
            ]
            sumstats: [route_meta, left_sumstats, right_sumstats]
            reference_ld_scores: [[id: meta.reference_bundle_id], reference_ld_scores]
            regression_weights: [[id: meta.reference_bundle_id], regression_weights]
        }

    LDSC_RG_LIABILITY(
        ch_ldsc_rg_liability.sumstats,
        ch_ldsc_rg_liability.reference_ld_scores,
        ch_ldsc_rg_liability.regression_weights,
    )

    def ch_ldsc_h2_native = LDSC_H2_OBSERVED.out.log
        .map { meta, observed_log -> [meta.request_id, meta, observed_log] }
        .join(
            LDSC_H2_LIABILITY.out.log.map { meta, liability_log -> [meta.request_id, liability_log] },
            remainder: true,
        )
        .join(ch_ldsc_h2_observed.munging_log, failOnDuplicate: true, failOnMismatch: true)
        .map { _request_id, meta, observed_log, liability_log, munging_log -> [meta, observed_log, liability_log ?: [], [munging_log]] }

    def ch_ldsc_rg_native = LDSC_RG_OBSERVED.out.log
        .map { meta, observed_log -> [meta.request_id, meta, observed_log] }
        .join(
            LDSC_RG_LIABILITY.out.log.map { meta, liability_log -> [meta.request_id, liability_log] },
            remainder: true,
        )
        .join(ch_ldsc_rg_observed.munging_logs, failOnDuplicate: true, failOnMismatch: true)
        .map { _request_id, meta, observed_log, liability_log, munging_logs -> [meta, observed_log, liability_log ?: [], munging_logs] }

    NORMALISE_LDSC(ch_ldsc_h2_native.mix(ch_ldsc_rg_native))

    //
    // SUBWORKFLOW: GCTA GREML heritability
    //
    // GCTA rejects a header row, so this route takes the headerless serialisations rather than the headered
    // ones the association routes use, and the trait sits at a fixed third column, which makes `--mpheno`
    // the constant 1 (set in conf/modules/gcta.config).
    //
    // The dense and LDMS matrix families retain distinct reuse keys and are adapted into the one public
    // GCTA heritability contract here. The middle GRM element is absent for GREML and is the MGRM manifest
    // for GREML-LDMS; the estimator selector makes the subworkflow enforce that distinction.
    def ch_greml_matrices = PREPARE_RELATEDNESS_MATRICES.out.gcta_dense
        .filter { meta, _grm_files -> !meta.relationship_id }
        .map { meta, grm_files -> [meta, [], grm_files, 'greml'] }
        .mix(
            PREPARE_RELATEDNESS_MATRICES.out.gcta_ldms
                .filter { meta, _mgrm, _grm_files -> !meta.relationship_id }
                .map { meta, mgrm, grm_files -> [meta, mgrm, grm_files, 'greml_ldms'] }
        )

    def ch_greml_inputs = ch_greml_matrices
        .combine(ch_gcta_phenotypes, by: 0)
        .multiMap { meta, mgrm, grm_files, estimator, phenotype, quant_covariates, cat_covariates ->
            def route_meta = meta + [gcta_estimator: estimator]
            grm: [route_meta, mgrm, grm_files]
            pheno: [route_meta, phenotype]
            qcovar: [route_meta, quant_covariates]
            covar: [route_meta, cat_covariates]
            estimator: [route_meta, estimator]
        }

    GRM_HERITABILITY_GCTA(
        ch_greml_inputs.grm,
        ch_greml_inputs.pheno,
        ch_greml_inputs.qcovar,
        ch_greml_inputs.covar,
        ch_greml_inputs.estimator,
    )

    //
    // PIPELINE ROUTE: primary dense GCTA bivariate REML relationship request
    //
    // The installed atomic component correctly requires the primary metadata ID to be the staged GRM
    // basename. Keep that native basename separate from request attribution and from the content-derived
    // matrix reuse key; all three identities reach the normalized provenance adapter.
    def ch_bivariate_matrices = PREPARE_RELATEDNESS_MATRICES.out.gcta_dense
        .filter { meta, _grm_files -> meta.relationship_id && meta.method == 'gcta_bivariate_reml' }
        .map { meta, grm_files ->
            def grm_id = grm_files.find { grm_file -> grm_file.name.endsWith('.grm.id') }
            if (!grm_id) {
                error("[nf-core/gwas] ERROR: pair request '${meta.request_id}' received a dense GCTA matrix without a .grm.id member")
            }
            def basename = grm_id.name.substring(0, grm_id.name.length() - '.grm.id'.length())
            [meta.request_id, meta + [matrix_basename: basename], grm_files]
        }

    def ch_prepared_relationships = PREPARE_BIVARIATE_TRAITS.out.phenotype
        .join(PREPARE_BIVARIATE_TRAITS.out.quant_covariates, remainder: true)
        .join(PREPARE_BIVARIATE_TRAITS.out.cat_covariates, remainder: true)
        .join(PREPARE_BIVARIATE_TRAITS.out.log, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, phenotype, quant_covariates, cat_covariates, pair_log ->
            [meta.relationship_id, phenotype, quant_covariates ?: [], cat_covariates ?: [], pair_log]
        }

    def ch_prepared_pairs = ch_relationships
        .map { meta, _genotype_files, _pair_quant_covariates, _pair_cat_covariates -> [meta.relationship_id, meta] }
        .combine(ch_prepared_relationships, by: 0)
        .map { _relationship_id, meta, phenotype, quant_covariates, cat_covariates, pair_log ->
            [meta.request_id, meta, phenotype, quant_covariates, cat_covariates, pair_log]
        }

    def ch_dense_prepared_pairs = ch_prepared_pairs
        .filter { _request_id, pair_meta, _phenotype, _quant_covariates, _cat_covariates, _pair_log -> pair_meta.method == 'gcta_bivariate_reml' }

    def ch_bivariate_invocations = ch_bivariate_matrices
        .join(ch_dense_prepared_pairs, failOnDuplicate: true, failOnMismatch: true)
        .multiMap { _request_id, matrix_meta, grm_files, pair_meta, phenotype, quant_covariates, cat_covariates, pair_log ->
            if (matrix_meta.relationship_id != pair_meta.relationship_id) {
                error("[nf-core/gwas] ERROR: pair request '${pair_meta.request_id}' matrix attribution disagrees with the prepared phenotype")
            }
            def route_meta = pair_meta + [
                id: matrix_meta.matrix_basename,
                matrix_key: matrix_meta.matrix_key,
                matrix_basename: matrix_meta.matrix_basename,
            ]
            grm: [route_meta, grm_files]
            pheno: [route_meta, phenotype, 1, 2]
            qcovar: [route_meta, quant_covariates]
            covar: [route_meta, cat_covariates]
            pair_log: [pair_meta.request_id, pair_log]
        }

    GCTA_BIVARIATEREML(
        ch_bivariate_invocations.grm,
        ch_bivariate_invocations.pheno,
        ch_bivariate_invocations.qcovar,
        ch_bivariate_invocations.covar,
    )

    def ch_dense_bivariate_native_results = GCTA_BIVARIATEREML.out.bivariate_results
        .map { meta, hsq -> [meta.request_id, meta, hsq] }
        .join(
            GCTA_BIVARIATEREML.out.log_file.map { meta, gcta_log -> [meta.request_id, gcta_log] },
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(ch_bivariate_invocations.pair_log, failOnDuplicate: true, failOnMismatch: true)
        .map { _request_id, meta, hsq, gcta_log, pair_log -> [meta, hsq, gcta_log, pair_log] }

    //
    // PIPELINE ROUTE: primary GCTA bivariate REML-LDMS relationship request
    //
    // The MGRM manifest basename is the installed atom's native identity. The request ID and the
    // matrix content key remain separate attribution fields so a unary GREML-LDMS request and a pair
    // request can share one scientifically identical matrix family without sharing result identity.
    def ch_bivariate_ldms_matrices = PREPARE_RELATEDNESS_MATRICES.out.gcta_ldms
        .filter { meta, _mgrm, _grm_files -> meta.relationship_id && 'gcta_bivariate_reml_ldms' == meta.method }
        .map { meta, mgrm, grm_files ->
            [meta.request_id, meta + [matrix_basename: mgrm.baseName], mgrm, grm_files]
        }

    def ch_ldms_prepared_pairs = ch_prepared_pairs
        .filter { _request_id, pair_meta, _phenotype, _quant_covariates, _cat_covariates, _pair_log -> pair_meta.method == 'gcta_bivariate_reml_ldms' }

    def ch_bivariate_ldms_invocations = ch_bivariate_ldms_matrices
        .join(ch_ldms_prepared_pairs, failOnDuplicate: true, failOnMismatch: true)
        .multiMap { _request_id, matrix_meta, mgrm, grm_files, pair_meta, phenotype, quant_covariates, cat_covariates, pair_log ->
            if (matrix_meta.relationship_id != pair_meta.relationship_id) {
                error("[nf-core/gwas] ERROR: pair request '${pair_meta.request_id}' LDMS matrix attribution disagrees with the prepared phenotype")
            }
            def route_meta = pair_meta + [
                id: matrix_meta.matrix_basename,
                matrix_key: matrix_meta.matrix_key,
                matrix_basename: matrix_meta.matrix_basename,
            ]
            mgrm: [route_meta, mgrm, grm_files]
            pheno: [route_meta, phenotype, 1, 2]
            qcovar: [route_meta, quant_covariates]
            covar: [route_meta, cat_covariates]
            pair_log: [pair_meta.request_id, pair_log]
        }

    GCTA_BIVARIATEREMLLDMS(
        ch_bivariate_ldms_invocations.mgrm,
        ch_bivariate_ldms_invocations.pheno,
        ch_bivariate_ldms_invocations.qcovar,
        ch_bivariate_ldms_invocations.covar,
    )

    def ch_ldms_bivariate_native_results = GCTA_BIVARIATEREMLLDMS.out.bivariate_results
        .map { meta, hsq -> [meta.request_id, meta, hsq] }
        .join(
            GCTA_BIVARIATEREMLLDMS.out.log_file.map { meta, gcta_log -> [meta.request_id, gcta_log] },
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(ch_bivariate_ldms_invocations.pair_log, failOnDuplicate: true, failOnMismatch: true)
        .map { _request_id, meta, hsq, gcta_log, pair_log -> [meta, hsq, gcta_log, pair_log] }

    def ch_bivariate_native_results = ch_dense_bivariate_native_results.mix(ch_ldms_bivariate_native_results)

    NORMALISE_GCTA_BIVARIATE(ch_bivariate_native_results)

    //
    // SUBWORKFLOWS: LDAK REML, Haseman-Elston and PCGC heritability
    //
    // Matrix construction and the per-analysis unrelated-subset routing are owned above by
    // PREPARE_RELATEDNESS_MATRICES. The three aliases preserve the reusable subworkflow's one-estimator
    // contract while allowing one analysis unit to select all three methods without changing its identity.
    // HE and PCGC additionally receive the numerical design built specifically for LDAK matrix adjustment;
    // the estimators themselves retain the original quantitative/categorical split.
    def ch_ldak_inputs = PREPARE_RELATEDNESS_MATRICES.out.ldak_kinship
        .join(ch_gcta_phenotypes, failOnDuplicate: true)
        .join(NORMALISE_PHENOTYPES.out.adjustment_covariates, remainder: true)
        .filter { record -> record.size() == 7 && record[1] != null }
        .map { meta, grm_files, keep, phenotype, quant_covariates, cat_covariates, adjustment_covariates ->
            [meta, grm_files, keep, phenotype, quant_covariates, cat_covariates, adjustment_covariates ?: []]
        }

    def ch_ldak_reml_inputs = ch_ldak_inputs
        .filter { meta, _grm_files, _keep, _phenotype, _quant_covariates, _cat_covariates, _adjustment_covariates ->
            'ldak_reml' in meta.heritability_methods
        }
        .multiMap { meta, grm_files, keep, phenotype, quant_covariates, cat_covariates, adjustment_covariates ->
            grm: [meta, grm_files]
            pheno: [meta, phenotype, meta.population_prevalence != null ? meta.population_prevalence : []]
            qcovar: [meta, quant_covariates]
            covar: [meta, cat_covariates]
            keep: [meta, keep ?: []]
            estimator: [meta, 'reml']
            adjustment_covar: [meta, adjustment_covariates]
        }

    GRM_HERITABILITY_LDAK_REML(
        ch_ldak_reml_inputs.grm,
        ch_ldak_reml_inputs.pheno,
        ch_ldak_reml_inputs.qcovar,
        ch_ldak_reml_inputs.covar,
        ch_ldak_reml_inputs.keep,
        ch_ldak_reml_inputs.estimator,
        ch_ldak_reml_inputs.adjustment_covar,
    )

    def ch_ldak_he_inputs = ch_ldak_inputs
        .filter { meta, _grm_files, _keep, _phenotype, _quant_covariates, _cat_covariates, _adjustment_covariates ->
            'ldak_he' in meta.heritability_methods
        }
        .multiMap { meta, grm_files, keep, phenotype, quant_covariates, cat_covariates, adjustment_covariates ->
            grm: [meta, grm_files]
            pheno: [meta, phenotype, []]
            qcovar: [meta, quant_covariates]
            covar: [meta, cat_covariates]
            keep: [meta, keep ?: []]
            estimator: [meta, 'he']
            adjustment_covar: [meta, adjustment_covariates]
        }

    GRM_HERITABILITY_LDAK_HE(
        ch_ldak_he_inputs.grm,
        ch_ldak_he_inputs.pheno,
        ch_ldak_he_inputs.qcovar,
        ch_ldak_he_inputs.covar,
        ch_ldak_he_inputs.keep,
        ch_ldak_he_inputs.estimator,
        ch_ldak_he_inputs.adjustment_covar,
    )

    def ch_ldak_pcgc_inputs = ch_ldak_inputs
        .filter { meta, _grm_files, _keep, _phenotype, _quant_covariates, _cat_covariates, _adjustment_covariates ->
            'ldak_pcgc' in meta.heritability_methods
        }
        .multiMap { meta, grm_files, keep, phenotype, quant_covariates, cat_covariates, adjustment_covariates ->
            grm: [meta, grm_files]
            pheno: [meta, phenotype, meta.population_prevalence]
            qcovar: [meta, quant_covariates]
            covar: [meta, cat_covariates]
            keep: [meta, keep ?: []]
            estimator: [meta, 'pcgc']
            adjustment_covar: [meta, adjustment_covariates]
        }

    // All constituent local modules report directly to the run-wide `versions` topic.
    GRM_HERITABILITY_LDAK_PCGC(
        ch_ldak_pcgc_inputs.grm,
        ch_ldak_pcgc_inputs.pheno,
        ch_ldak_pcgc_inputs.qcovar,
        ch_ldak_pcgc_inputs.covar,
        ch_ldak_pcgc_inputs.keep,
        ch_ldak_pcgc_inputs.estimator,
        ch_ldak_pcgc_inputs.adjustment_covar,
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
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def multiqc_analysis_plan = file("${projectDir}/assets/multiqc_analysis_plan.yml", checkIfExists: true)
    def ch_analysis_plan = ch_analysis_metadata.map { analysis_metadata -> analysisPlanJson(multiqc_analysis_plan, analysis_metadata) }
    ch_multiqc_files = ch_multiqc_files.mix(ch_analysis_plan.collectFile(name: 'analysis_plan_mqc.json', sort: true))
    def summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = ch_method_metadata.map { method_metadata ->
        def analysis_metadata = method_metadata.findAll { record -> record.domain == 'analysis' }.collect { record -> record.meta }
        def summary_unary_metadata = method_metadata.findAll { record -> record.domain == 'summary_unary' }.collect { record -> record.meta }
        def relationship_metadata = method_metadata.findAll { record -> record.domain == 'pairwise' }.collect { record -> record.meta }
        def selected_methods = [
            association: analysis_metadata.collectMany { meta -> meta.association_methods }.unique().sort(),
            heritability: (analysis_metadata.collectMany { meta -> meta.heritability_methods } + summary_unary_metadata.collect { meta -> meta.method }).unique().sort(),
            pairwise: relationship_metadata.collectMany { meta -> meta.relationship_methods }.unique().sort(),
        ]
        methodsDescriptionText(multiqc_custom_methods_description, selected_methods)
    }
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'gwas'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo
                    ? file(multiqc_logo, checkIfExists: true)
                    : file("${projectDir}/assets/nf-core-gwas_logo_light.png", checkIfExists: true),
                [],
                [],
            ]
        }
    )

    emit:
    canonical_summary_statistics  = CANONICALISE_SUMMARY_STATISTICS.out.summary_statistics // channel: [ val(meta), path(canonical_summary_statistics) ]
    summary_statistics_provenance = CANONICALISE_SUMMARY_STATISTICS.out.provenance // channel: [ val(meta), path(provenance) ]
    multiqc_report                = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: [ [ path(report) ] ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def getLdakSummaryArguments(meta) {
    def args = new ArrayList(meta.native_args ?: [])
    def option_names = args.findAll { token -> token instanceof String && token.startsWith('--') }.collect { token -> token.split('=', 2)[0] }
    if (!option_names.contains('--cutoff') && !option_names.contains('--truncate')) {
        args.addAll(['--cutoff', '0.01'])
    }
    if (meta.method == 'ldak_sumher' && meta.is_binary && meta.population_prevalence != null && meta.sample_prevalence != null) {
        args.addAll(['--prevalence', meta.population_prevalence.toString(), '--ascertainment', meta.sample_prevalence.toString()])
    }
    if (meta.method == 'ldak_sumcors' && meta.left_is_binary && meta.right_is_binary && meta.left_population_prevalence != null && meta.left_sample_prevalence != null && meta.right_population_prevalence != null && meta.right_sample_prevalence != null) {
        args.addAll(
            [
                '--prevalence',
                meta.left_population_prevalence.toString(),
                '--ascertainment',
                meta.left_sample_prevalence.toString(),
                '--prevalence2',
                meta.right_population_prevalence.toString(),
                '--ascertainment2',
                meta.right_sample_prevalence.toString(),
            ]
        )
    }
    return args
}

def getLdakSummaryRuntime() {
    return 'ghcr.io/lyh970817/gwas/ldak:6.3-b755ab7@sha256:f2b2157559e4346cab5e9f478ab70fc76359743ef06522fed9ad23769d735a6e'
}
