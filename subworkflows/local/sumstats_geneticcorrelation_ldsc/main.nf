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
    // Each endpoint travels as a named record so the fit reads its parts by key rather than by position.
    def ch_endpoints = ch_pairs
        .flatMap { meta, meta2, reference_ld_scores, meta3, regression_weights, left_producer_id, right_producer_id ->
            [[left_producer_id, 0], [right_producer_id, 1]].collect { producer_id, ordinal ->
                [
                    producer_id,
                    [
                        meta: meta,
                        reference_meta: meta2,
                        reference: reference_ld_scores,
                        weights_meta: meta3,
                        weights: regression_weights,
                        ordinal: ordinal,
                    ],
                ]
            }
        }
        .groupTuple()
        .join(ch_munged, failOnMismatch: true, failOnDuplicate: true)
        .flatMap { _producer_id, endpoints, sumstats ->
            endpoints.collect { endpoint -> [groupKey(endpoint.meta.id, 2), endpoint + [sumstats: sumstats]] }
        }
    // Ordinal 0 is the endpoint the pair declared first, so it carries the pair identity and both resources.
    def ch_fits = ch_endpoints
        .groupTuple()
        .map { _pair_key, endpoints -> endpoints.sort { endpoint -> endpoint.ordinal } }
        .multiMap { left, right ->
            sumstats: [left.meta, left.sumstats, right.sumstats]
            reference: [left.reference_meta, left.reference]
            weights: [left.weights_meta, left.weights]
        }
    LDSC_RG(ch_fits.sumstats, ch_fits.reference, ch_fits.weights)

    emit:
    munged_sumstats = LDSC_MUNGESUMSTATS.out.munged_sumstats // channel: [ val(producer_meta), path(sumstats.gz) ]
    munging_log     = LDSC_MUNGESUMSTATS.out.log // channel: [ val(producer_meta), path(log) ]
    rg_log          = LDSC_RG.out.log // channel: [ val(pair_meta), path(log) ]
}
