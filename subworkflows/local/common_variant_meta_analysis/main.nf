// Run explicitly selected common-variant meta-analysis models on one genome-wide request.

include { GWASLAB_METAANALYZE       } from '../../../modules/local/gwaslab/metaanalyze/main'
include { METASOFT_RE2              } from '../../../modules/local/metasoft/re2/main'
include { CUSTOM_PREPAREMRMEGAINPUT } from '../../../modules/local/custom/preparemrmegainput/main'
include { MRMEGA                    } from '../../../modules/local/mrmega/main'

workflow COMMON_VARIANT_META_ANALYSIS {
    take:
    ch_request // channel: [ val(meta), path(parents), val(study_names), val(input_format), val(genome_build), val(trait_type), val(models) ]

    main:
    GWASLAB_METAANALYZE(
        ch_request.map { meta, parents, study_names, input_format, genome_build, _trait_type, models ->
            [meta, parents, study_names, input_format, genome_build, models.contains('random')]
        }
    )

    ch_selection = ch_request.map { meta, _parents, study_names, _input_format, _genome_build, trait_type, models ->
        [meta, study_names, trait_type, models]
    }

    ch_metasoft = GWASLAB_METAANALYZE.out.metasoft_input
        .join(ch_selection, failOnDuplicate: true, failOnMismatch: true)
        .filter { _meta, _matrix, _study_names, _trait_type, models -> models.contains('re2') }
        .map { meta, matrix, _study_names, _trait_type, _models -> [meta, matrix] }
    METASOFT_RE2(ch_metasoft)

    ch_mrmega = GWASLAB_METAANALYZE.out.mrmega_input
        .join(ch_selection, failOnDuplicate: true, failOnMismatch: true)
        .filter { _meta, _aligned, _study_names, _trait_type, models -> models.contains('mrmega') }

    CUSTOM_PREPAREMRMEGAINPUT(
        ch_mrmega.map { meta, aligned, study_names, trait_type, _models ->
            [meta, aligned, study_names, trait_type]
        }
    )

    // One trait_type selects both the exported effect columns and MR-MEGA's native input-column mode.
    MRMEGA(
        CUSTOM_PREPAREMRMEGAINPUT.out.study_files.join(CUSTOM_PREPAREMRMEGAINPUT.out.filelist, failOnDuplicate: true, failOnMismatch: true).join(
            ch_mrmega.map { meta, _aligned, _study_names, trait_type, _models -> [meta, trait_type] },
            failOnDuplicate: true,
            failOnMismatch: true,
        )
    )

    emit:
    fixed           = GWASLAB_METAANALYZE.out.fixed // channel: [ val(meta), path(fixed) ]
    gwaslab_log     = GWASLAB_METAANALYZE.out.log // channel: [ val(meta), path(log) ]
    random_effects  = GWASLAB_METAANALYZE.out.random_effects // channel: [ val(meta), path(random) ]
    metasoft_result = METASOFT_RE2.out.result // channel: [ val(meta), path(result) ]
    metasoft_log    = METASOFT_RE2.out.log // channel: [ val(meta), path(log) ]
    mrmega_result   = MRMEGA.out.result // channel: [ val(meta), path(result) ]
    mrmega_log      = MRMEGA.out.log // channel: [ val(meta), path(log) ]
}
