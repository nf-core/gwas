// Run explicitly selected common-variant meta-analysis models on one genome-wide request.

include { GWASLAB_META_ANALYZE } from '../../../modules/local/gwaslab/meta_analyze/main'
include { METASOFT_RE2 } from '../../../modules/local/metasoft/re2/main'
include { PREPARE_MRMEGA_INPUT } from '../../../modules/local/prepare_mrmega_input/main'
include { MRMEGA } from '../../../modules/local/mrmega/main'

workflow COMMON_VARIANT_META_ANALYSIS {
    take:
    ch_request // channel: [ val(meta), path(parents), val(study_names), val(input_format), val(genome_build), val(trait_type), val(models), val(axes) ]

    main:
    GWASLAB_META_ANALYZE(
        ch_request.map { meta, parents, study_names, input_format, genome_build, _trait_type, models, _axes ->
            [meta, parents, study_names, input_format, genome_build, models.contains('random')]
        }
    )

    ch_selection = ch_request.map { meta, _parents, study_names, _input_format, _genome_build, trait_type, models, axes ->
        [meta, study_names, trait_type, models, axes]
    }

    ch_metasoft = GWASLAB_META_ANALYZE.out.metasoft_input
        .join(ch_selection)
        .filter { _meta, _matrix, _study_names, _trait_type, models, _axes -> models.contains('re2') }
        .map { meta, matrix, _study_names, _trait_type, _models, _axes -> [meta, matrix] }
    METASOFT_RE2(ch_metasoft)

    ch_mrmega = GWASLAB_META_ANALYZE.out.mrmega_input
        .join(ch_selection)
        .filter { _meta, _aligned, _study_names, _trait_type, models, _axes -> models.contains('mrmega') }

    PREPARE_MRMEGA_INPUT(
        ch_mrmega.map { meta, aligned, study_names, trait_type, _models, _axes ->
            [meta, aligned, study_names, trait_type]
        }
    )

    ch_mrmega_context = ch_mrmega.map { meta, _aligned, _study_names, trait_type, _models, axes ->
        [meta, axes, trait_type]
    }

    MRMEGA(
        PREPARE_MRMEGA_INPUT.out.study_files.join(PREPARE_MRMEGA_INPUT.out.filelist).join(ch_mrmega_context).map { meta, study_files, filelist, axes, trait_type ->
            [meta, study_files, filelist, axes, trait_type]
        }
    )

    emit:
    fixed = GWASLAB_META_ANALYZE.out.fixed // channel: [ val(meta), path(fixed) ]
    gwaslab_log = GWASLAB_META_ANALYZE.out.log // channel: [ val(meta), path(log) ]
    random_effects = GWASLAB_META_ANALYZE.out.random_effects // channel: [ val(meta), path(random) ]
    metasoft_result = METASOFT_RE2.out.result // channel: [ val(meta), path(result) ]
    metasoft_log = METASOFT_RE2.out.log // channel: [ val(meta), path(log) ]
    mrmega_result = MRMEGA.out.result // channel: [ val(meta), path(result) ]
    mrmega_log = MRMEGA.out.log // channel: [ val(meta), path(log) ]
}
