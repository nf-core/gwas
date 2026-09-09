// Build one LD-by-MAF component plan: the LD scores of a focal PLINK 1 bundle and the ordered, disjoint SNP
// groups derived from them. The plan is the shared input of every LD- and MAF-stratified matrix family, so it
// is built here once rather than inside any one tool's matrix builder.
// Every constituent process reports directly to the run-wide versions topic, so this subworkflow emits no versions.

// MODULE: Local to the pipeline
include { GCTA_CALCULATELDSCORES      } from '../../../modules/local/gcta/calculateldscores/main'
include { CUSTOM_GCTASTRATIFYLDSCORES } from '../../../modules/local/custom/gctastratifyldscores/main'

workflow PLINK_PREPARE_LDMS_PLAN_GCTA {
    take:
    ch_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam) ], focal genotypes
    ch_ld_score_region_kb // channel: [ val(meta2), val(ld_score_region_kb) ], LD-score region size
    ch_ld_bins // channel: [ val(meta3), val(ld_bins) ], LD-score bin count
    ch_maf_edges // channel: [ val(meta4), val(maf_edges) ], MAF edges

    main:
    ch_plan_inputs = ch_genotypes
        .map { meta, bed, bim, fam -> tuple(meta.id, tuple(meta, bed, bim, fam)) }
        .join(ch_ld_score_region_kb.map { meta2, ld_score_region_kb -> tuple(meta2.id, ld_score_region_kb) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_ld_bins.map { meta3, ld_bins -> tuple(meta3.id, ld_bins) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_maf_edges.map { meta4, maf_edges -> tuple(meta4.id, maf_edges) }, by: 0, failOnDuplicate: true, failOnMismatch: true)

    ch_ld_score_inputs = ch_plan_inputs.multiMap { _plan_id, genotypes, ld_score_region_kb, _ld_bins, _maf_edges ->
        genotypes: genotypes
        region_kb: ld_score_region_kb
    }
    GCTA_CALCULATELDSCORES(ch_ld_score_inputs.genotypes, ch_ld_score_inputs.region_kb)

    ch_stratify_inputs = GCTA_CALCULATELDSCORES.out.ld_scores
        .map { meta, ld_scores -> tuple(meta.id, meta, ld_scores) }
        .join(ch_plan_inputs.map { plan_id, _genotypes, _ld_score_region_kb, ld_bins, maf_edges -> tuple(plan_id, ld_bins, maf_edges) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _plan_id, meta, ld_scores, ld_bins, maf_edges -> tuple(meta, ld_scores, ld_bins, maf_edges) }
    CUSTOM_GCTASTRATIFYLDSCORES(ch_stratify_inputs)

    emit:
    plan      = CUSTOM_GCTASTRATIFYLDSCORES.out.strata_bundle // channel: [ val(meta), path(strata_manifest), [ path(snp_group_file), ... ] ], LD-major MAF-minor row order
    ld_scores = GCTA_CALCULATELDSCORES.out.ld_scores // channel: [ val(meta), path(ld_scores) ]
}
