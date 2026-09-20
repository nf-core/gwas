// Prepare declared summary producers once and fit each explicitly referenced heritability request.
include { LDSC_MUNGESUMSTATS } from '../../../modules/local/ldsc/mungesumstats/main'
include { LDSC_H2            } from '../../../modules/local/ldsc/h2/main'

workflow SUMSTATS_HERITABILITY_LDSC {
    take:
    ch_producers // channel: [ val(meta), path(sumstats), val(meta2), path(merge_alleles) ]; optional merge_alleles: []
    ch_requests // channel: [ val(meta), val(meta2), path(reference_ld_scores), val(reference_prefix), val(meta3), path(regression_weights), val(weights_prefix), val(producer_id) ]

    main:
    def ch_munge = ch_producers.multiMap { meta, sumstats, meta2, merge_alleles ->
        sumstats: [meta, sumstats]
        alleles: [meta2, merge_alleles]
    }
    LDSC_MUNGESUMSTATS(ch_munge.sumstats, ch_munge.alleles)

    def ch_munged = LDSC_MUNGESUMSTATS.out.munged_sumstats.map { meta, sumstats -> [meta.id, sumstats] }
    def ch_fits = ch_requests
        .map { meta, meta2, reference_ld_scores, reference_prefix, meta3, regression_weights, weights_prefix, producer_id ->
            [producer_id, [meta, meta2, reference_ld_scores, reference_prefix, meta3, regression_weights, weights_prefix]]
        }
        .groupTuple()
        .join(ch_munged, failOnMismatch: true, failOnDuplicate: true)
        .flatMap { _producer_id, requests, sumstats -> requests.collect { request -> request + [sumstats] } }
        .multiMap { meta, meta2, reference_ld_scores, reference_prefix, meta3, regression_weights, weights_prefix, sumstats ->
            sumstats: [meta, sumstats]
            reference: [meta2, reference_ld_scores, reference_prefix]
            weights: [meta3, regression_weights, weights_prefix]
        }
    LDSC_H2(ch_fits.sumstats, ch_fits.reference, ch_fits.weights)

    emit:
    munged_sumstats = LDSC_MUNGESUMSTATS.out.munged_sumstats // channel: [ val(producer_meta), path(sumstats.gz) ]
    munging_log     = LDSC_MUNGESUMSTATS.out.log // channel: [ val(producer_meta), path(log) ]
    h2_log          = LDSC_H2.out.log // channel: [ val(request_meta), path(log) ]
}
