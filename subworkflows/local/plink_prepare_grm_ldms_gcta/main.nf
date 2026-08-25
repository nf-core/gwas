// Build LD-by-MAF-stratified GCTA genetic relationship matrices
// Every constituent process reports directly to the run-wide versions topic, so this subworkflow emits no versions.

// MODULE: Local to the pipeline
include { GCTA_CALCULATELDSCORES        } from '../../../modules/local/gcta/calculateldscores/main'
include { CUSTOM_GCTASTRATIFYLDSCORES   } from '../../../modules/local/custom/gctastratifyldscores/main'
include { CUSTOM_GCTACREATEMGRMMANIFEST } from '../../../modules/local/custom/gctacreatemgrmmanifest/main'

// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
include { PLINK_PREPARE_GRM_GCTA        } from '../plink_prepare_grm_gcta/main'

workflow PLINK_PREPARE_GRM_LDMS_GCTA {
    take:
    ch_genotypes // channel: [ val(meta), path(mfile), path(bed), path(bim), path(fam) ], focal genotypes
    ch_ld_score_region_kb // channel: [ val(meta2), val(ld_score_region_kb) ], LD-score region size
    ch_ld_bins // channel: [ val(meta3), val(ld_bins) ], LD-score bin count
    ch_maf_edges // channel: [ val(meta4), val(maf_edges) ], MAF edges
    ch_n_parts // channel: [ val(meta5), val(requested_parts) ], part count

    main:
    ch_analysis_inputs = ch_genotypes
        .map { meta, mfile, bed, bim, fam -> tuple(meta.id, tuple(meta, mfile, bed, bim, fam)) }
        .join(ch_ld_score_region_kb.map { meta2, ld_score_region_kb -> tuple(meta2.id, ld_score_region_kb) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_ld_bins.map { meta3, ld_bins -> tuple(meta3.id, ld_bins) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_maf_edges.map { meta4, maf_edges -> tuple(meta4.id, maf_edges) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_n_parts.map { meta5, requested_parts -> tuple(meta5.id, requested_parts) }, by: 0, failOnDuplicate: true, failOnMismatch: true)

    ch_ld_score_inputs = ch_analysis_inputs.multiMap { _analysis_id, genotypes, ld_score_region_kb, _ld_bins, _maf_edges, _requested_parts ->
        genotypes: tuple(genotypes[0], genotypes[2], genotypes[3], genotypes[4])
        region_kb: ld_score_region_kb
    }
    GCTA_CALCULATELDSCORES(ch_ld_score_inputs.genotypes, ch_ld_score_inputs.region_kb)

    ch_stratify_inputs = GCTA_CALCULATELDSCORES.out.ld_scores
        .map { meta, ld_scores -> tuple(meta.id, meta, ld_scores) }
        .join(ch_analysis_inputs.map { analysis_id, _genotypes, _ld_score_region_kb, ld_bins, maf_edges, _requested_parts -> tuple(analysis_id, ld_bins, maf_edges) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _analysis_id, meta, ld_scores, ld_bins, maf_edges -> tuple(meta, ld_scores, ld_bins, maf_edges) }
    CUSTOM_GCTASTRATIFYLDSCORES(ch_stratify_inputs)

    ch_stratum_state = CUSTOM_GCTASTRATIFYLDSCORES.out.strata_bundle
        .map { meta, manifest, snp_group_files -> tuple(meta.id, meta, manifest, snp_group_files) }
        .join(ch_analysis_inputs.map { analysis_id, genotypes, _ld_score_region_kb, _ld_bins, _maf_edges, requested_parts -> tuple(analysis_id, genotypes, requested_parts) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .flatMap { _analysis_id, focal_meta, manifest, snp_group_files, genotypes, requested_parts ->
            def rows = manifest.readLines()
            def columns = rows[0].split('\t', -1).toList()
            def groups = (snp_group_files instanceof List ? snp_group_files : [snp_group_files]).collectEntries { snp_group_file -> [snp_group_file.name, snp_group_file] }
            rows
                .drop(1)
                .withIndex()
                .collect { row, index ->
                    def provenance = [columns, row.split('\t', -1).toList()].transpose().collectEntries()
                    def ordinal = index + 1
                    tuple([id: "${focal_meta.id}__${String.format('%06d', ordinal)}"], focal_meta, ordinal, rows.size() - 1, provenance, genotypes, groups[provenance.group_filename], requested_parts)
                }
        }

    ch_dense_inputs = ch_stratum_state.multiMap { work_meta, focal_meta, _ordinal, stratum_count, provenance, genotypes, snp_group_file, requested_parts ->
        genotypes: tuple(work_meta, genotypes[1], genotypes[2], genotypes[3], genotypes[4])
        snp_group: tuple(work_meta, snp_group_file)
        n_parts: tuple(work_meta, requested_parts)
        restore: tuple(work_meta.id, focal_meta, stratum_count)
        strata: tuple(focal_meta, provenance.model_key, provenance.stratum_key, provenance)
    }
    PLINK_PREPARE_GRM_GCTA(ch_dense_inputs.genotypes, ch_dense_inputs.snp_group, ch_dense_inputs.n_parts)

    ch_mgrm_inputs = PLINK_PREPARE_GRM_GCTA.out.grm_files
        .map { work_meta, grm_files -> tuple(work_meta.id, grm_files) }
        .join(ch_dense_inputs.restore, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _work_id, grm_files, focal_meta, stratum_count -> tuple(groupKey(focal_meta, stratum_count), grm_files) }
        .groupTuple()
        .map { key, grm_file_lists -> tuple(key.getGroupTarget(), grm_file_lists.flatten()) }
    CUSTOM_GCTACREATEMGRMMANIFEST(ch_mgrm_inputs)

    emit:
    mgrm_bundle = CUSTOM_GCTACREATEMGRMMANIFEST.out.mgrm_bundle // channel: [ val(meta), path(mgrm), path(grm_files) ]
    strata      = ch_dense_inputs.strata // channel: [ val(meta), val(model_key), val(stratum_key), val(provenance) ], manifest order
}
