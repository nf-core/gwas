// Build one LDAK kinship matrix and optionally derive its unrelated-sample subset.
// Constituent local modules report versions directly to the run-wide topic, so this subworkflow emits no versions.
//

// MODULE: Local to the pipeline
include { LDAK_CALCKINS } from '../../../modules/local/ldak/calckins/main'
include { LDAK_FILTER   } from '../../../modules/local/ldak/filter/main'
include { LDAK_SUBGRM   } from '../../../modules/local/ldak/subgrm/main'

workflow PLINK_PREPARE_GRM_LDAK {
    take:
    ch_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam), val(power) ], once per analysis
    ch_weights // channel: [ val(meta2), path(weights_file) ], use [] for the optional file, once per analysis
    ch_filter_relatedness // channel: [ val(meta3), val(filter_relatedness) ], route selector once per analysis

    main:
    def ch_calckins_inputs = ch_genotypes
        .map { meta, bed, bim, fam, power -> tuple(meta.id, tuple(meta, bed, bim, fam, power)) }
        .join(
            ch_weights.map { meta2, weights_file -> tuple(meta2.id, tuple(meta2, weights_file)) },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .join(
            ch_filter_relatedness.map { meta3, filter_relatedness -> tuple(meta3.id, filter_relatedness) },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .multiMap { analysis_id, genotypes, weights, filter_relatedness ->
            genotypes: genotypes
            weights: weights
            routing: tuple(analysis_id, genotypes[0], filter_relatedness)
        }

    LDAK_CALCKINS(ch_calckins_inputs.genotypes, ch_calckins_inputs.weights)

    def ch_unfiltered_grm = LDAK_CALCKINS.out.ldak_grm.map { meta, grm_bin, grm_id, grm_details, grm_adjust ->
        tuple(meta, [grm_bin, grm_id, grm_details, grm_adjust])
    }
    def ch_grm_routes = ch_unfiltered_grm
        .map { meta, grm_files -> tuple(meta.id, tuple(meta, grm_files)) }
        .join(
            ch_calckins_inputs.routing,
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .branch { _analysis_id, grm, focal_meta, filter_relatedness ->
            direct: !filter_relatedness
            return tuple(focal_meta, grm[1])
            filtered: filter_relatedness
            return tuple(focal_meta, grm[1])
        }

    LDAK_FILTER(ch_grm_routes.filtered)

    def ch_subgrm_state = ch_grm_routes.filtered
        .map { meta, grm_files -> tuple(meta.id, tuple(meta, grm_files)) }
        .join(
            LDAK_FILTER.out.filtered_list.map { meta, keep, _lose -> tuple(meta.id, keep) },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .multiMap { analysis_id, grm, keep ->
            subgrm: tuple([id: "${analysis_id}.unrelated"], grm[1], keep)
            focal_meta: tuple("${analysis_id}.unrelated", grm[0])
        }
    LDAK_SUBGRM(ch_subgrm_state.subgrm)

    def ch_filtered_analysis_grm = LDAK_SUBGRM.out.sub_grm
        .map { execution_meta, grm_bin, grm_id, grm_details, grm_adjust ->
            tuple(execution_meta.id, [grm_bin, grm_id, grm_details, grm_adjust])
        }
        .join(
            ch_subgrm_state.focal_meta,
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { _execution_id, grm_files, focal_meta -> tuple(focal_meta, grm_files) }
    def ch_analysis_grm = ch_grm_routes.direct.mix(ch_filtered_analysis_grm)

    emit:
    analysis_grm   = ch_analysis_grm // channel: [ val(meta), path(grm_files) ], selected collected GRM bundle
    unfiltered_grm = ch_unfiltered_grm // channel: [ val(meta), path(grm_files) ], collected all-sample GRM bundle
    filtered_list  = LDAK_FILTER.out.filtered_list // channel: [ val(meta), path(keep), path(lose) ], filtered-route records only
    maxrel         = LDAK_FILTER.out.maxrel // channel: [ val(meta), path(maxrel) ], optional filtered-route records only
}
