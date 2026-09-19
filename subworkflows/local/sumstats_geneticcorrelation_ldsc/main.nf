// Munge each declared endpoint once and fit ordered pairs that reference those producers.
include { LDSC_MUNGESUMSTATS } from '../../../modules/local/ldsc/mungesumstats/main'
include { LDSC_RG            } from '../../../modules/local/ldsc/rg/main'

workflow SUMSTATS_GENETICCORRELATION_LDSC {
    take:
    ch_producers // channel: [ val(meta), path(sumstats), val(meta2), path(merge_alleles) ]; optional merge_alleles: []
    ch_pairs // channel: [ val(meta), val(meta2), path(reference_ld_scores), val(meta3), path(regression_weights), val(left_producer_id), val(right_producer_id) ]

    main:
    def ch_munge = ch_producers.multiMap { meta, sumstats, meta2, merge_alleles ->
        sumstats: [meta, sumstats]
        alleles: [meta2, merge_alleles]
    }
    LDSC_MUNGESUMSTATS(ch_munge.sumstats, ch_munge.alleles)

    def ch_munged = LDSC_MUNGESUMSTATS.out.munged_sumstats.map { meta, sumstats -> [meta.id, sumstats] }
    // Resolve both roles together so each producer joins once even when it is shared by several pairs.
    def ch_endpoints = ch_pairs
        .flatMap { meta, meta2, reference_ld_scores, meta3, regression_weights, left_producer_id, right_producer_id ->
            [[left_producer_id, 0], [right_producer_id, 1]].collect { producer_id, ordinal ->
                [producer_id, [meta, meta2, reference_ld_scores, meta3, regression_weights, ordinal]]
            }
        }
        .groupTuple()
        .join(ch_munged, failOnMismatch: true, failOnDuplicate: true)
        .flatMap { _producer_id, endpoints, sumstats ->
            endpoints.collect { meta, meta2, reference_ld_scores, meta3, regression_weights, ordinal ->
                [groupKey(meta.id, 2), [meta, meta2, reference_ld_scores, meta3, regression_weights, ordinal, sumstats]]
            }
        }
    def ch_fits = ch_endpoints
        .groupTuple()
        .map { _pair_key, endpoints -> endpoints.sort { endpoint -> endpoint[5] } }
        .multiMap { left, right ->
            sumstats: [left[0], left[6], right[6]]
            reference: [left[1], left[2]]
            weights: [left[3], left[4]]
        }
    LDSC_RG(ch_fits.sumstats, ch_fits.reference, ch_fits.weights)

    emit:
    munged_sumstats = LDSC_MUNGESUMSTATS.out.munged_sumstats // channel: [ val(producer_meta), path(sumstats.gz) ]
    munging_log     = LDSC_MUNGESUMSTATS.out.log // channel: [ val(producer_meta), path(log) ]
    rg_log          = LDSC_RG.out.log // channel: [ val(pair_meta), path(log) ]
}
