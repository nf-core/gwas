// Construct and merge dense GCTA genetic relationship matrix parts.
// Both processes report on the run-wide versions topic, so this subworkflow emits no versions.

// MODULE: Local to the pipeline
include { GCTA_MAKEGRMPART         } from '../../../modules/local/gcta/makegrmpart/main'
include { CUSTOM_GCTAMERGEGRMPARTS } from '../../../modules/local/custom/gctamergegrmparts/main'

workflow PLINK_PREPARE_GRM_GCTA {
    take:
    ch_genotypes // channel: [ val(meta), path(mfile), path(bed_pgen), path(bim_pvar), path(fam_psam) ], mandatory
    ch_snp_group_file // channel: [ val(meta2), path(snp_group_file) ], optional; use [] when absent
    ch_n_parts // channel: [ val(meta3), val(requested_parts) ], mandatory

    main:
    ch_make_jobs = ch_genotypes
        .map { meta, mfile, bed_pgen, bim_pvar, fam_psam -> [meta.id, meta, mfile, bed_pgen, bim_pvar, fam_psam] }
        .join(
            ch_snp_group_file.map { meta2, snp_group_file -> [meta2.id, meta2, snp_group_file] },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_n_parts.map { meta3, requested_parts -> [meta3.id, requested_parts] },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .flatMap { analysis_id, meta, mfile, bed_pgen, bim_pvar, fam_psam, meta2, snp_group_file, requested_parts ->
            def part_width = requested_parts.toString().size()
            (1..requested_parts).collect { part_index ->
                [analysis_id, meta, requested_parts, String.format("%0${part_width}d", part_index), mfile, bed_pgen, bim_pvar, fam_psam, meta2, snp_group_file]
            }
        }
    ch_make_inputs = ch_make_jobs.multiMap { _analysis_id, meta, requested_parts, padded_part, mfile, bed_pgen, bim_pvar, fam_psam, meta2, snp_group_file ->
        genotypes: [meta, requested_parts, padded_part, mfile, bed_pgen, bim_pvar, fam_psam]
        snp_group: [meta2, snp_group_file]
    }
    GCTA_MAKEGRMPART(ch_make_inputs.genotypes, ch_make_inputs.snp_group)

    ch_gathered_parts = GCTA_MAKEGRMPART.out.grm_files
        .map { meta, grm_part_files, nparts_gcta, _part_gcta_job ->
            [groupKey(meta, nparts_gcta), grm_part_files]
        }
        .groupTuple()
        .map { key, grm_part_file_lists -> [key.getGroupTarget(), grm_part_file_lists.flatten(), grm_part_file_lists.size()] }

    CUSTOM_GCTAMERGEGRMPARTS(ch_gathered_parts.map { meta, grm_part_files, _effective_parts -> [meta, grm_part_files] })

    emit:
    grm_files       = CUSTOM_GCTAMERGEGRMPARTS.out.grm_files // channel: [ val(meta), path(grm_files) ]
    effective_parts = ch_gathered_parts.map { meta, _grm_part_files, effective_parts -> [meta, effective_parts] } // channel: [ val(meta), val(effective_parts) ]
}
