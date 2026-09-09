// Build LD-by-MAF-stratified GCTA genetic relationship matrices from a completed component plan.
// The plan itself is built once by `plink_prepare_ldms_plan_gcta` and passed in, because a second tool's
// stratified matrix family consumes the same ordered manifest and the same SNP group files.
// Every constituent process reports directly to the run-wide versions topic, so this subworkflow emits no versions.

// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
include { PLINK_PREPARE_GRM_GCTA } from '../plink_prepare_grm_gcta/main'

workflow PLINK_PREPARE_GRM_LDMS_GCTA {
    take:
    ch_plan // channel: [ val(meta), path(strata_manifest), [ path(snp_group_file), ... ] ], ordered component plan
    ch_genotypes // channel: [ val(meta2), path(mfile), path(bed), path(bim), path(fam) ], focal genotypes
    ch_n_parts // channel: [ val(meta3), val(requested_parts) ], part count

    main:
    ch_analysis_inputs = ch_plan
        .map { meta, strata_manifest, snp_group_files -> tuple(meta.id, meta, strata_manifest, snp_group_files) }
        .join(ch_genotypes.map { meta2, mfile, bed, bim, fam -> tuple(meta2.id, tuple(meta2, mfile, bed, bim, fam)) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_n_parts.map { meta3, requested_parts -> tuple(meta3.id, requested_parts) }, by: 0, failOnDuplicate: true, failOnMismatch: true)

    ch_stratum_state = ch_analysis_inputs.flatMap { _analysis_id, focal_meta, manifest, snp_group_files, genotypes, requested_parts ->
        readLdmsComponents(manifest, snp_group_files).collect { component ->
            tuple(
                [id: "${focal_meta.id}__${String.format('%06d', component.ordinal)}"],
                focal_meta,
                component.ordinal,
                component.count,
                component.provenance,
                genotypes,
                component.group_file,
                requested_parts,
            )
        }
    }

    ch_dense_inputs = ch_stratum_state.multiMap { work_meta, focal_meta, ordinal, stratum_count, provenance, genotypes, snp_group_file, requested_parts ->
        genotypes: tuple(work_meta, genotypes[1], genotypes[2], genotypes[3], genotypes[4])
        snp_group: tuple(work_meta, snp_group_file)
        n_parts: tuple(work_meta, requested_parts)
        restore: tuple(work_meta.id, focal_meta, ordinal, stratum_count)
        strata: tuple(focal_meta, provenance.model_key, provenance.stratum_key, provenance)
    }
    PLINK_PREPARE_GRM_GCTA(ch_dense_inputs.genotypes, ch_dense_inputs.snp_group, ch_dense_inputs.n_parts)

    ch_grm_families = PLINK_PREPARE_GRM_GCTA.out.grm_files
        .map { work_meta, grm_files -> tuple(work_meta.id, grm_files) }
        .join(ch_dense_inputs.restore, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _work_id, grm_files, focal_meta, ordinal, stratum_count -> tuple(groupKey(focal_meta, stratum_count), ordinal, grm_files) }
        .groupTuple()
        .map { key, ordinals, grm_file_lists ->
            def ordered_strata = orderLdmsStrata(ordinals, grm_file_lists)
            def grm_prefixes = ordered_strata.collect { _ordinal, grm_files ->
                def grm_id_name = grm_files.find { grm_file -> grm_file.name.endsWith('.grm.id') }.name
                grm_id_name.substring(0, grm_id_name.length() - '.grm.id'.length())
            }
            tuple(key.getGroupTarget(), ordered_strata.collectMany { _ordinal, grm_files -> grm_files }, grm_prefixes)
        }

    emit:
    grm_family = ch_grm_families // channel: [ val(meta), path(grm_files), val(grm_prefixes) ], explicit stratum order
    strata     = ch_dense_inputs.strata // channel: [ val(meta), val(model_key), val(stratum_key), val(provenance) ], declared non-empty stratum order
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def orderLdmsStrata(ordinals, grm_file_lists) {
    return [ordinals, grm_file_lists]
        .transpose()
        .sort { left, right -> left[0] <=> right[0] }
}

// Expand a completed component plan into ordered component records. The ordinal is the 1-based manifest row,
// and the manifest row order is LD-major and MAF-minor with empty strata omitted, which is the stratifier's
// contract. `count` travels with every record so the gather above can size its `groupKey` without re-reading
// the file.
//
// Deliberately duplicated rather than imported from the plan builder: each composition that expands a plan
// carries its own copy, so no component depends on a Groovy function living inside another submission
// candidate. Every copy is covered by its own `tests/main.function.nf.test`.
def readLdmsComponents(manifest, snp_group_files) {
    def rows = manifest.readLines()
    def columns = rows[0].split('\t', -1).toList()
    def groups = (snp_group_files instanceof List ? snp_group_files : [snp_group_files]).collectEntries { snp_group_file -> [snp_group_file.name, snp_group_file] }
    return rows
        .drop(1)
        .withIndex()
        .collect { row, index ->
            def provenance = [columns, row.split('\t', -1).toList()].transpose().collectEntries()
            [
                ordinal: index + 1,
                count: rows.size() - 1,
                stratum_key: provenance.stratum_key,
                predictor_count: provenance.predictor_count as int,
                group_file: groups[provenance.group_filename],
                provenance: provenance,
            ]
        }
}
