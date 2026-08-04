// Prepare each distinct cohort once as canonical PLINK 2 and derive PLINK 1 only for routes that need it.
// All three processes report on the run-wide versions topic, so this subworkflow emits no versions.

// MODULE: Local to the pipeline
include { PLINK2_MAKEPGEN } from '../../../modules/local/plink2/makepgen/main'
include { PLINK2_MAKEBED  } from '../../../modules/local/plink2/makebed/main'
include { PLINK2_VCF      } from '../../../modules/local/plink2/vcf/main'

workflow PREPARE_COHORT_GENOTYPES {
    take:
    ch_analyses // channel: [ val(meta), [ path(genotype_file), ... ] ]

    main:

    //
    // The prepared-key set is a reduction of the analysis-key set, not a second collection compared
    // against it. Every cohort key emitted here came from an analysis unit, so the prepared set is
    // total; `unique` collapses the repeats, so it is unique. Both properties hold structurally and
    // there is nothing left to assert afterwards — which matters, because the `combine` below
    // reports neither a missing key nor a duplicate one. Building the cohort set from any other
    // source, or asserting totality after the fact, would reintroduce exactly the failure this
    // shape removes.
    //
    // Collapsing on the cohort key alone is safe because input validation has already refused any
    // samplesheet whose rows share a `cohort_id` while naming different genotype files, so every
    // discarded repeat is known to be identical to the one kept.
    //
    // The cohort key is carried only so `unique` can collapse on it, so it is dropped again here,
    // in one place, the moment it has served that purpose. [ cohort_meta, [ genotype_file, ... ] ]
    def ch_cohorts = ch_analyses
        .map { meta, genotype_files ->
            [meta.cohort, [id: meta.cohort, genotype_format: meta.genotype_format], genotype_files]
        }
        .unique { cohort_id, _cohort_meta, _genotype_files -> cohort_id }
        .map { _cohort_id, cohort_meta, genotype_files -> [cohort_meta, genotype_files] }

    def ch_by_format = ch_cohorts.branch { cohort_meta, _genotype_files ->
        plink2: cohort_meta.genotype_format == 'plink2'
        plink1: cohort_meta.genotype_format == 'plink1'
        vcf: cohort_meta.genotype_format == 'vcf'
    }

    //
    // PLINK 2 is the canonical representation, so a cohort supplied in it is already prepared and is
    // passed through untouched: converting it would cost a round trip and produce the same bundle.
    // Nothing the pipeline built is therefore published for such a cohort under the save control —
    // the canonical bundle is the researcher's own input, already on disk where they put it.
    //
    def ch_supplied = ch_by_format.plink2.map { cohort_meta, genotype_files ->
        def (pgen, psam, pvar) = genotype_files
        [cohort_meta, pgen, psam, pvar]
    }

    //
    // MODULE: Convert a PLINK 1 cohort to the canonical PLINK 2 bundle
    //
    PLINK2_MAKEPGEN(
        ch_by_format.plink1.map { cohort_meta, genotype_files ->
            def (bed, bim, fam) = genotype_files
            [cohort_meta, bed, bim, fam]
        }
    )

    //
    // MODULE: Convert a VCF cohort to the canonical PLINK 2 bundle
    //
    PLINK2_VCF(
        ch_by_format.vcf.map { cohort_meta, genotype_files -> [cohort_meta, genotype_files.first()] }
    )

    def ch_cohort_genotypes = ch_supplied
    ch_cohort_genotypes = ch_cohort_genotypes.mix(PLINK2_MAKEPGEN.out.pgen)
    ch_cohort_genotypes = ch_cohort_genotypes.mix(PLINK2_VCF.out.pgen)

    //
    // The cohort-to-analysis seam. `cohort_id` to `analysis_id` is one-to-many, so this is a
    // combining operator rather than a join: `join` emits one pair per key and so silently drops
    // every analysis unit after the first on each cohort, while its `failOnDuplicate` form raises on
    // the second instead. `combine(by: 0)` pairs every analysis unit with its cohort's bundle, and
    // because exactly one bundle exists per cohort key the product is exactly one element per
    // analysis unit, carrying the focal analysis meta rather than the cohort meta.
    //
    def ch_genotypes = ch_analyses
        .map { meta, _genotype_files -> [meta.cohort, meta] }
        .combine(
            ch_cohort_genotypes.map { cohort_meta, pgen, psam, pvar -> [cohort_meta.id, pgen, psam, pvar] },
            by: 0
        )
        .map { _cohort_id, meta, pgen, psam, pvar -> [meta, pgen, psam, pvar] }

    //
    // LDAK and GCTA GREML-LDMS consume PLINK 1. The derivative is lazy — only a cohort with a route
    // that needs it is converted — and cohort-keyed, so all such analyses share one conversion. It is
    // always derived from the canonical PLINK 2 bundle, even when the researcher supplied PLINK 1, so
    // every input encoding crosses the same compatibility seam.
    //
    def ch_plink1_analyses = ch_analyses.filter { meta, _genotype_files ->
        'ldak_kvik' in (meta.association_methods ?: []) || (meta.heritability_methods ?: []).any { method -> method in ['ldak_reml', 'ldak_he', 'ldak_pcgc'] } || 'gcta_greml_ldms' in (meta.heritability_methods ?: [])
    }

    def ch_plink1_cohort_genotypes = ch_plink1_analyses
        .map { meta, _genotype_files -> [meta.cohort] }
        .unique { cohort_id -> cohort_id }
        .combine(
            ch_cohort_genotypes.map { cohort_meta, pgen, psam, pvar -> [cohort_meta.id, cohort_meta, pgen, psam, pvar] },
            by: 0
        )
        .map { _cohort_id, cohort_meta, pgen, psam, pvar -> [cohort_meta, pgen, psam, pvar] }

    PLINK2_MAKEBED(ch_plink1_cohort_genotypes)

    // The cohort-to-analysis relationship is one-to-many here for the same reason as the canonical
    // fan-out above. Replace the conversion's cohort metadata with the unchanged focal analysis meta.
    def ch_plink1_genotypes = ch_plink1_analyses
        .map { meta, _genotype_files -> [meta.cohort, meta] }
        .combine(
            PLINK2_MAKEBED.out.bed.map { cohort_meta, bed, bim, fam -> [cohort_meta.id, bed, bim, fam] },
            by: 0
        )
        .map { _cohort_id, meta, bed, bim, fam -> [meta, bed, bim, fam] }

    emit:
    genotypes        = ch_genotypes // channel: [ val(meta), path(pgen), path(psam), path(pvar) ]
    cohort_genotypes = ch_cohort_genotypes // channel: [ val(cohort_meta), path(pgen), path(psam), path(pvar) ]
    plink1_genotypes = ch_plink1_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam) ], analyses needing PLINK 1
}
