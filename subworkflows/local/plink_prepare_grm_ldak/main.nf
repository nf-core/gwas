// Build one LDAK base kinship and optionally derive its unrelated-sample subset.
// Constituent local modules report versions directly to the run-wide topic, so this subworkflow emits no versions.

// MODULE: Local to the pipeline
include { LDAK_CALCKINS } from '../../../modules/local/ldak/calckins/main'
include { LDAK_FILTER   } from '../../../modules/local/ldak/filter/main'
include { LDAK_SUBGRM   } from '../../../modules/local/ldak/subgrm/main'

workflow PLINK_PREPARE_GRM_LDAK {
    take:
    ch_genotypes // channel: [ val(base_meta), path(bed), path(bim), path(fam), val(power) ]
    ch_weights // channel: [ val(base_meta2), path(weights_file) ], use [] for the optional file
    ch_subset_requests // channel: [ val(parent_key), val(derived_meta) ], only requested subset derivatives

    main:
    def ch_calckins_inputs = ch_genotypes
        .map { meta, bed, bim, fam, power -> [meta.id, [meta, bed, bim, fam, power]] }
        .join(
            ch_weights.map { meta2, weights_file -> [meta2.id, [meta2, weights_file]] },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .multiMap { _base_id, genotypes, weights ->
            genotypes: genotypes
            weights: weights
        }

    LDAK_CALCKINS(ch_calckins_inputs.genotypes, ch_calckins_inputs.weights)

    def ch_base_grm = LDAK_CALCKINS.out.ldak_grm.map { meta, grm_bin, grm_id, grm_details, grm_adjust ->
        [meta, [grm_bin, grm_id, grm_details, grm_adjust]]
    }

    // A child request selects its completed parent by key and replaces the base execution metadata with the
    // derivative metadata before either native transformation runs.
    def ch_subset_inputs = ch_base_grm
        .map { base_meta, grm_files -> [base_meta.key, grm_files] }
        .join(ch_subset_requests, by: 0)
        .map { _parent_key, grm_files, derived_meta -> [derived_meta, grm_files] }

    LDAK_FILTER(ch_subset_inputs)

    def ch_subgrm_inputs = ch_subset_inputs
        .map { derived_meta, grm_files -> [derived_meta.id, derived_meta, grm_files] }
        .join(
            LDAK_FILTER.out.filtered_list.map { derived_meta, keep, _lose -> [derived_meta.id, keep] },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { _derived_id, derived_meta, grm_files, keep -> [derived_meta, grm_files, keep] }

    LDAK_SUBGRM(ch_subgrm_inputs)

    def ch_subset_grm = LDAK_SUBGRM.out.sub_grm.map { derived_meta, grm_bin, grm_id, grm_details, grm_adjust ->
        [derived_meta, [grm_bin, grm_id, grm_details, grm_adjust]]
    }

    emit:
    base_grm      = ch_base_grm // channel: [ val(base_meta), path(grm_files) ]
    subset_grm    = ch_subset_grm // channel: [ val(derived_meta), path(grm_files) ]
    filtered_list = LDAK_FILTER.out.filtered_list // channel: [ val(derived_meta), path(keep), path(lose) ]
    maxrel        = LDAK_FILTER.out.maxrel // channel: [ val(derived_meta), path(maxrel) ], optional
}
