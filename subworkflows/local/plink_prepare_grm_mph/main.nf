// Build one native MPH genetic relationship matrix over the whole declared variant universe.
// The universe is expressed as a weight column rather than as a genotype subset: MPH treats a variant weighted
// 0 exactly as it treats an omitted one, so the SNP-information file stays auditable against the BIM line for
// line while still excluding what the caller declared out of scope.
// Both processes report directly to the run-wide versions topic, so this subworkflow emits no versions.

// MODULE: Local to the pipeline
include { CUSTOM_MPHSNPINFO } from '../../../modules/local/custom/mphsnpinfo/main'
include { MPH_MAKEGRM       } from '../../../modules/local/mph/makegrm/main'

workflow PLINK_PREPARE_GRM_MPH {
    take:
    ch_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam) ], focal genotypes
    autosome_count // value: highest integer chromosome code counted as autosomal

    main:
    CUSTOM_MPHSNPINFO(
        ch_genotypes.map { meta, _bed, bim, _fam -> tuple(meta, bim) },
        ch_genotypes.map { meta, _bed, _bim, _fam -> tuple(meta, [], [], ['all']) },
        autosome_count,
    )

    ch_make_inputs = ch_genotypes
        .map { meta, bed, bim, fam -> tuple(meta.id, meta, bed, bim, fam) }
        .join(
            CUSTOM_MPHSNPINFO.out.snp_info.map { meta, snp_info, weight_names -> tuple(meta.id, snp_info, weight_names) },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .multiMap { _analysis_id, meta, bed, bim, fam, snp_info, weight_names ->
            genotypes: tuple(meta, bed, bim, fam)
            // The weight column comes from the adapter's echoed list rather than from a literal here, so the
            // one-component name is stated in exactly one place.
            snp_info: tuple(meta, snp_info, weight_names.first())
        }

    MPH_MAKEGRM(ch_make_inputs.genotypes, ch_make_inputs.snp_info)

    emit:
    grm_files = MPH_MAKEGRM.out.grm_files // channel: [ val(meta), path(grm_files) ], the [.grm.bin, .grm.iid] bundle
    counts    = CUSTOM_MPHSNPINFO.out.counts // channel: [ val(meta), path(counts_tsv) ]
}
