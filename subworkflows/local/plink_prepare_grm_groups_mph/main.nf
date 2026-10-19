// Portable dependencies: named SNP-group conversion and native MPH matrix construction.
include { CUSTOM_MPHSNPINFO } from '../../../modules/local/custom/mphsnpinfo/main'
include { MPH_MAKEGRM       } from '../../../modules/local/mph/makegrm/main'

workflow PLINK_PREPARE_GRM_GROUPS_MPH {
    take:
    ch_groups // channel: [ val(meta), path(group_manifest), [ path(snp_group_file), ... ] ], ordered named groups
    ch_genotypes // channel: [ val(meta2), path(bed), path(bim), path(fam) ], focal genotypes

    main:
    ch_family_inputs = ch_groups
        .map { meta, group_manifest, snp_group_files -> tuple(meta.id, meta, group_manifest, snp_group_files) }
        .join(ch_genotypes.map { meta2, bed, bim, fam -> tuple(meta2.id, bed, bim, fam) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { analysis_id, meta, group_manifest, snp_group_files, bed, bim, fam ->
            tuple(analysis_id, meta, group_manifest, snp_group_files, bed, bim, fam, readSnpGroups(group_manifest))
        }

    ch_snp_info_inputs = ch_family_inputs.multiMap { _analysis_id, meta, group_manifest, snp_group_files, _bed, bim, _fam, components ->
        bim: tuple(meta, bim)
        plan: tuple(meta, group_manifest, snp_group_files, components.collect { component -> component.group_id })
    }
    CUSTOM_MPHSNPINFO(ch_snp_info_inputs.bim, ch_snp_info_inputs.plan)

    ch_component_state = ch_family_inputs
        .map { analysis_id, meta, _group_manifest, _snp_group_files, bed, bim, fam, components -> tuple(analysis_id, meta, bed, bim, fam, components) }
        .join(CUSTOM_MPHSNPINFO.out.snp_info.map { meta, snp_info, weight_names -> tuple(meta.id, snp_info, weight_names) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .flatMap { _analysis_id, focal_meta, bed, bim, fam, components, snp_info, weight_names ->
            components.collect { component ->
                // MPH uses matrix prefixes as variance-component labels, so retain the group ID.
                tuple(
                    [id: "${focal_meta.id}__${String.format('%06d', component.ordinal)}_${component.group_id}"],
                    focal_meta,
                    component.ordinal,
                    component.count,
                    bed,
                    bim,
                    fam,
                    snp_info,
                    weight_names[component.ordinal - 1],
                )
            }
        }

    ch_make_inputs = ch_component_state.multiMap { work_meta, focal_meta, ordinal, component_count, bed, bim, fam, snp_info, weight_name ->
        genotypes: tuple(work_meta, bed, bim, fam)
        snp_info: tuple(focal_meta, snp_info, weight_name)
        restore: tuple(work_meta.id, focal_meta, ordinal, component_count)
    }
    MPH_MAKEGRM(ch_make_inputs.genotypes, ch_make_inputs.snp_info)

    ch_grm_families = MPH_MAKEGRM.out.grm_files
        .map { work_meta, grm_files -> tuple(work_meta.id, grm_files) }
        .join(ch_make_inputs.restore, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _work_id, grm_files, focal_meta, ordinal, component_count -> tuple(groupKey(focal_meta, component_count), ordinal, grm_files) }
        .groupTuple()
        .map { key, ordinals, grm_file_lists ->
            def ordered_components = orderSnpGroups(ordinals, grm_file_lists)
            def grm_prefixes = ordered_components.collect { _ordinal, grm_files ->
                def grm_iid_name = grm_files.find { grm_file -> grm_file.name.endsWith('.grm.iid') }.name
                grm_iid_name.substring(0, grm_iid_name.length() - '.grm.iid'.length())
            }
            tuple(key.getGroupTarget(), ordered_components.collectMany { _ordinal, grm_files -> grm_files }, grm_prefixes)
        }

    emit:
    grm_family = ch_grm_families // channel: [ val(meta), path(grm_files), val(grm_prefixes) ], explicit component order
    counts     = CUSTOM_MPHSNPINFO.out.counts // channel: [ val(meta), path(counts_tsv) ], one per family
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def orderSnpGroups(ordinals, grm_file_lists) {
    return [ordinals, grm_file_lists]
        .transpose()
        .sort { left, right -> left[0] <=> right[0] }
}

// Row order defines group order; the caller owns the scientific grouping and selected SNP universe.
def readSnpGroups(manifest) {
    def rows = manifest.readLines()
    def columns = rows[0].split('\t', -1).toList()
    return rows
        .drop(1)
        .withIndex()
        .collect { row, index ->
            def group = [columns, row.split('\t', -1).toList()].transpose().collectEntries()
            [
                ordinal: index + 1,
                count: rows.size() - 1,
                group_id: group.group_id,
            ]
        }
}
