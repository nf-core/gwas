// Construct and merge dense GCTA genetic relationship matrix parts.
// Both processes report on the run-wide versions topic, so this subworkflow emits no versions.

// MODULE: Local to the pipeline
include { GCTA_MAKEGRMPART   } from '../../../modules/local/gcta/makegrmpart/main'
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
            (1..requested_parts).collect { part_index ->
                [analysis_id, meta, requested_parts, part_index, mfile, bed_pgen, bim_pvar, fam_psam, meta2, snp_group_file]
            }
        }
    ch_make_inputs = ch_make_jobs.multiMap { _analysis_id, meta, requested_parts, part_index, mfile, bed_pgen, bim_pvar, fam_psam, meta2, snp_group_file ->
        genotypes: [meta, requested_parts, part_index, mfile, bed_pgen, bim_pvar, fam_psam]
        snp_group: [meta2, snp_group_file]
    }
    GCTA_MAKEGRMPART(ch_make_inputs.genotypes, ch_make_inputs.snp_group)

    ch_gathered_parts = GCTA_MAKEGRMPART.out.grm_files
        .map { meta, grm_part_files, nparts_gcta, part_gcta_job ->
            [groupKey(meta, nparts_gcta), part_gcta_job as int, grm_part_files]
        }
        .groupTuple()
        .map { key, part_indices, grm_part_file_lists ->
            def ordered_parts = orderGrmPartRecords(part_indices, grm_part_file_lists)
            def ordered_files = ordered_parts.collectMany { part -> [part.bin, part.n_bin, part.id] }
            def staged_records = ordered_parts.collect { part -> [part: part.part, bin: part.bin.name, n_bin: part.n_bin.name, id: part.id.name] }
            [key.getGroupTarget(), ordered_files, staged_records]
        }

    CUSTOM_GCTAMERGEGRMPARTS(ch_gathered_parts)

    emit:
    grm_files = CUSTOM_GCTAMERGEGRMPARTS.out.grm_files // channel: [ val(meta), path(grm_files) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def orderGrmPartRecords(part_indices, grm_part_file_lists) {
    return [part_indices, grm_part_file_lists]
        .transpose()
        .sort { left, right -> left[0] <=> right[0] }
        .collect { part_index, grm_part_files ->
            [
                part: part_index,
                bin: grm_part_files.find { grm_part_file -> grm_part_file.name.endsWith('.grm.bin') },
                n_bin: grm_part_files.find { grm_part_file -> grm_part_file.name.endsWith('.grm.N.bin') },
                id: grm_part_files.find { grm_part_file -> grm_part_file.name.endsWith('.grm.id') },
            ]
        }
}
