// Export one MR-MEGA study file per aligned summary-statistics source, then run one genome-wide MR-MEGA meta-regression.
include { CUSTOM_PREPAREMRMEGAINPUT } from '../../../modules/local/custom/preparemrmegainput/main'
include { MRMEGA                    } from '../../../modules/local/mrmega/main'

workflow SUMSTATS_METAANALYSIS_MRMEGA {
    take:
    ch_aligned // channel: [ val(meta), path(aligned_studies), val(study_names), val(trait_type) ]

    main:
    CUSTOM_PREPAREMRMEGAINPUT(ch_aligned)

    // One trait_type selects both the exported effect columns and MR-MEGA's native input-column mode.
    MRMEGA(
        CUSTOM_PREPAREMRMEGAINPUT.out.study_files
            .join(CUSTOM_PREPAREMRMEGAINPUT.out.filelist, failOnDuplicate: true, failOnMismatch: true)
            .join(ch_aligned.map { meta, _aligned, _study_names, trait_type -> [meta, trait_type] }, failOnDuplicate: true, failOnMismatch: true)
    )

    emit:
    result     = MRMEGA.out.result // channel: [ val(meta), path(result) ]
    mrmega_log = MRMEGA.out.log // channel: [ val(meta), path(log) ]
}
