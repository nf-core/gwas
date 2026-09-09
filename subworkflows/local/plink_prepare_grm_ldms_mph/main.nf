// Build an explicitly ordered LD-by-MAF-stratified family of native MPH genetic relationship matrices from a
// completed component plan. One SNP-information file carries every component as its own weight column, and one
// matrix is built per column from the full genotype bundle: passing the whole bundle with a weight column and
// extracting the component's variants first are byte-identical, so there is no per-component extract step.
// Every constituent process reports directly to the run-wide versions topic, so this subworkflow emits no versions.

// MODULE: Local to the pipeline
include { CUSTOM_MPHSNPINFO } from '../../../modules/local/custom/mphsnpinfo/main'
include { MPH_MAKEGRM       } from '../../../modules/local/mph/makegrm/main'

workflow PLINK_PREPARE_GRM_LDMS_MPH {
    take:
    ch_plan // channel: [ val(meta), path(strata_manifest), [ path(snp_group_file), ... ] ], ordered component plan
    ch_genotypes // channel: [ val(meta2), path(bed), path(bim), path(fam) ], focal genotypes
    autosome_count // value: highest integer chromosome code counted as autosomal

    main:
    ch_family_inputs = ch_plan
        .map { meta, strata_manifest, snp_group_files -> tuple(meta.id, meta, strata_manifest, snp_group_files) }
        .join(ch_genotypes.map { meta2, bed, bim, fam -> tuple(meta2.id, bed, bim, fam) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        // The plan is a completed process output here, so expanding it in an operator is legal; the same read
        // inside a `script:` block is not, because a staged `path` is a stage-name handle at that point.
        .map { analysis_id, meta, strata_manifest, snp_group_files, bed, bim, fam ->
            tuple(analysis_id, meta, strata_manifest, snp_group_files, bed, bim, fam, readLdmsComponents(strata_manifest, snp_group_files))
        }

    ch_snp_info_inputs = ch_family_inputs.multiMap { _analysis_id, meta, strata_manifest, snp_group_files, _bed, bim, _fam, components ->
        bim: tuple(meta, bim)
        plan: tuple(meta, strata_manifest, snp_group_files, components.collect { component -> component.stratum_key })
    }
    CUSTOM_MPHSNPINFO(ch_snp_info_inputs.bim, ch_snp_info_inputs.plan, autosome_count)

    ch_component_state = ch_family_inputs
        .map { analysis_id, meta, _strata_manifest, _snp_group_files, bed, bim, fam, components -> tuple(analysis_id, meta, bed, bim, fam, components) }
        .join(CUSTOM_MPHSNPINFO.out.snp_info.map { meta, snp_info, weight_names -> tuple(meta.id, snp_info, weight_names) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .flatMap { _analysis_id, focal_meta, bed, bim, fam, components, snp_info, weight_names ->
            components.collect { component ->
                // The stratum key rides in the work identifier as well as the ordinal because MPH copies the
                // `--grm_list` prefix verbatim into the result's `vc_name`. Without it the primary published
                // result would name its components by an opaque reuse key and the GCTA-to-MPH component
                // mapping would be purely positional. The `%06d` ordinal stays first so the gather below and
                // the GCTA-side naming convention are unchanged.
                tuple(
                    [id: "${focal_meta.id}__${String.format('%06d', component.ordinal)}_${component.stratum_key}"],
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
        // The SNP-information file is one artefact of the whole family rather than of this component, so it
        // carries the family's identity. That is also what keeps the atom's two-role tag readable instead of
        // repeating one long component identifier twice.
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
            def ordered_components = orderLdmsStrata(ordinals, grm_file_lists)
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
// Deliberately duplicated rather than imported from the plan builder or the GCTA builder: each composition
// that expands a plan carries its own copy, so no component depends on a Groovy function living inside another
// submission candidate. Every copy is covered by its own `tests/main.function.nf.test`.
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
