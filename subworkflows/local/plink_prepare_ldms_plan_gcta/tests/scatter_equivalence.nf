// Test harness for the per-chromosome scatter in PLINK_PREPARE_LDMS_PLAN_GCTA. Beside the composition it runs
// one genome-wide GCTA LD-score pass and the stratifier on that pass, on the same bundles and settings, so a
// test can compare what the two routes produce.

include { PLINK_PREPARE_LDMS_PLAN_GCTA                                        } from '../main'
include { GCTA_CALCULATELDSCORES as GCTA_CALCULATELDSCORES_GENOMEWIDE           } from '../../../../modules/local/gcta/calculateldscores/main'
include { CUSTOM_GCTASTRATIFYLDSCORES as CUSTOM_GCTASTRATIFYLDSCORES_GENOMEWIDE } from '../../../../modules/local/custom/gctastratifyldscores/main'

workflow LDMS_PLAN_SCATTER_EQUIVALENCE {
    take:
    ch_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam) ]
    ch_ld_score_region_kb // channel: [ val(meta2), val(ld_score_region_kb) ]
    ch_ld_bins // channel: [ val(meta3), val(ld_bins) ]
    ch_maf_edges // channel: [ val(meta4), val(maf_edges) ]

    main:
    PLINK_PREPARE_LDMS_PLAN_GCTA(ch_genotypes, ch_ld_score_region_kb, ch_ld_bins, ch_maf_edges)

    GCTA_CALCULATELDSCORES_GENOMEWIDE(
        ch_genotypes
            .map { meta, bed, bim, fam -> tuple(meta.id, meta, bed, bim, fam) }
            .join(ch_ld_score_region_kb.map { meta2, ld_score_region_kb -> tuple(meta2.id, ld_score_region_kb) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
            .map { _plan_id, meta, bed, bim, fam, ld_score_region_kb -> tuple(meta, bed, bim, fam, ld_score_region_kb, []) }
    )
    CUSTOM_GCTASTRATIFYLDSCORES_GENOMEWIDE(
        GCTA_CALCULATELDSCORES_GENOMEWIDE.out.ld_scores
            .map { meta, ld_scores -> tuple(meta.id, meta, ld_scores) }
            .join(ch_ld_bins.map { meta3, ld_bins -> tuple(meta3.id, ld_bins) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
            .join(ch_maf_edges.map { meta4, maf_edges -> tuple(meta4.id, maf_edges) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
            .map { _plan_id, meta, ld_scores, ld_bins, maf_edges -> tuple(meta, ld_scores, ld_bins, maf_edges) }
    )

    emit:
    scattered_ld_scores   = PLINK_PREPARE_LDMS_PLAN_GCTA.out.ld_scores // channel: [ val(meta), path(ld_scores) ]
    scattered_plan        = PLINK_PREPARE_LDMS_PLAN_GCTA.out.plan // channel: [ val(meta), path(strata_manifest), [ path(snp_group_file), ... ] ]
    genome_wide_ld_scores = GCTA_CALCULATELDSCORES_GENOMEWIDE.out.ld_scores // channel: [ val(meta), path(ld_scores) ]
    genome_wide_plan      = CUSTOM_GCTASTRATIFYLDSCORES_GENOMEWIDE.out.strata_bundle // channel: [ val(meta), path(strata_manifest), [ path(snp_group_file), ... ] ]
}
