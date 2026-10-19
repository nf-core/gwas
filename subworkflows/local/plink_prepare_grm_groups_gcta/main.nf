// Build an ordered GCTA matrix family from one or more caller-defined SNP groups.
// Row partitioning and ordered merging remain inside the independently callable single-matrix builder.

// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
include { PLINK_PREPARE_GRM_GCTA } from '../plink_prepare_grm_gcta/main'

workflow PLINK_PREPARE_GRM_GROUPS_GCTA {
    take:
    ch_groups // channel: [ val(meta), path(group_manifest), [ path(snp_group_file), ... ] ], ordered named-group table and referenced files
    ch_genotypes // channel: [ val(meta2), path(mfile), path(bed), path(bim), path(fam) ], focal genotypes; meta2.id equals meta.id
    ch_n_parts // channel: [ val(meta3), val(requested_parts) ], part count; meta3.id equals meta.id

    main:
    ch_analysis_inputs = ch_groups
        .map { meta, group_manifest, snp_group_files -> tuple(meta.id, meta, group_manifest, snp_group_files) }
        .join(ch_genotypes.map { meta2, mfile, bed, bim, fam -> tuple(meta2.id, tuple(meta2, mfile, bed, bim, fam)) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_n_parts.map { meta3, requested_parts -> tuple(meta3.id, requested_parts) }, by: 0, failOnDuplicate: true, failOnMismatch: true)

    // The declared group label is carried into the per-group work identity, so it survives into the native
    // matrix basenames and into the public `grm_prefixes` the caller uses to name variance components. The
    // leading ordinal keeps the key collision-free and sortable when two labels differ only in case.
    ch_group_state = ch_analysis_inputs.flatMap { _analysis_id, focal_meta, manifest, snp_group_files, genotypes, requested_parts ->
        readNamedGroups(manifest, snp_group_files).collect { component ->
            tuple(
                [id: "${focal_meta.id}__${String.format('%06d', component.ordinal)}_${component.group_id}"],
                focal_meta,
                component.ordinal,
                component.count,
                genotypes,
                component.group_file,
                requested_parts,
            )
        }
    }

    ch_dense_inputs = ch_group_state.multiMap { work_meta, focal_meta, ordinal, group_count, genotypes, snp_group_file, requested_parts ->
        genotypes: tuple(work_meta, genotypes[1], genotypes[2], genotypes[3], genotypes[4])
        snp_group: tuple(work_meta, snp_group_file)
        n_parts: tuple(work_meta, requested_parts)
        restore: tuple(work_meta.id, focal_meta, ordinal, group_count)
    }
    PLINK_PREPARE_GRM_GCTA(ch_dense_inputs.genotypes, ch_dense_inputs.snp_group, ch_dense_inputs.n_parts)

    ch_grm_families = PLINK_PREPARE_GRM_GCTA.out.grm_files
        .map { work_meta, grm_files -> tuple(work_meta.id, grm_files) }
        .join(ch_dense_inputs.restore, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _work_id, grm_files, focal_meta, ordinal, group_count -> tuple(groupKey(focal_meta, group_count), ordinal, grm_files) }
        .groupTuple()
        .map { key, ordinals, grm_file_lists ->
            def ordered_groups = orderGroupMatrices(ordinals, grm_file_lists)
            def grm_prefixes = ordered_groups.collect { _ordinal, grm_files ->
                def grm_id_name = grm_files.find { grm_file -> grm_file.name.endsWith('.grm.id') }.name
                grm_id_name.substring(0, grm_id_name.length() - '.grm.id'.length())
            }
            tuple(key.getGroupTarget(), ordered_groups.collectMany { _ordinal, grm_files -> grm_files }, grm_prefixes)
        }

    emit:
    grm_family = ch_grm_families // channel: [ val(meta), path(grm_files), val(grm_prefixes) ], declared group order
}

def orderGroupMatrices(ordinals, grm_file_lists) {
    return [ordinals, grm_file_lists]
        .transpose()
        .sort { left, right -> left[0] <=> right[0] }
}

// The table order defines matrix order independently of file arrival or task completion.
def readNamedGroups(manifest, snp_group_files) {
    def rows = manifest.readLines()
    def columns = rows[0].split('\t', -1).toList()
    def groups = (snp_group_files instanceof List ? snp_group_files : [snp_group_files]).collectEntries { group_file -> [group_file.name, group_file] }
    return rows
        .drop(1)
        .withIndex()
        .collect { row, index ->
            def record = [columns, row.split('\t', -1).toList()].transpose().collectEntries()
            [ordinal: index + 1, count: rows.size() - 1, group_id: record.group_id, group_file: groups[record.group_filename]]
        }
}
