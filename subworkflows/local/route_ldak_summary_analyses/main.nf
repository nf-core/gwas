// Route nf-core/gwas summary-scale LDAK requests — SumHer heritability and SumCors genetic correlation —
// from GWASLab-standard summary statistics through the required LDAK-format adapter to native results.
// SumHer and SumCors are one controller because they share that preparation and endpoint resolution.
// This is pipeline routing, reuse identity and scientific/publication policy, not an nf-core/modules
// submission candidate. Every constituent process reports directly to the run-wide versions topic, so this
// subworkflow emits no versions, and it reads no params, no workflow and no projectDir.

// MODULES: Upstream-ready components used inside a pipeline-local route
include { LDAK_SUMCORS                    } from '../../../modules/local/ldak/sumcors/main'
include { LDAK_SUMHER                     } from '../../../modules/local/ldak/sumher/main'

// MODULE: Local to the pipeline
include { PREPARE_LDAK_SUMMARY_STATISTICS } from '../../../modules/local/prepare_ldak_summary_statistics/main'

workflow ROUTE_LDAK_SUMMARY_ANALYSES {
    take:
    ch_sumher_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ], the ldak_sumher unary requests; LDAK-family bundles supply [] for elements 2-4
    ch_sumcors_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ], the ldak_sumcors pair requests, same convention
    ch_summary_statistics // channel: [ val(meta), path(gwaslab_summary_statistics) ], one element per distinct summary_statistics_id

    main:

    // Adapt each distinct summary once, regardless of how many unary, pairwise or named sensitivity
    // requests consume it. The adapter owns only the deterministic GWASLab-to-LDAK column transform;
    // each request retains its own tagging reference, effective native arguments and publication identity.
    def ch_ldak_requested_summary_ids = ch_sumher_requests
        .map { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file -> [meta.summary_statistics_id] }
        .mix(
            ch_sumcors_requests.flatMap { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file ->
                [[meta.left_summary_statistics_id], [meta.right_summary_statistics_id]]
            }
        )
        .unique()

    def ch_ldak_summaries = ch_ldak_requested_summary_ids
        .combine(
            ch_summary_statistics.map { meta, summary_statistics -> [meta.summary_statistics_id, meta, summary_statistics] },
            by: 0
        )
        .map { _summary_statistics_id, meta, summary_statistics -> [meta, summary_statistics] }

    PREPARE_LDAK_SUMMARY_STATISTICS(ch_ldak_summaries)

    def ch_prepared_ldak_summaries = PREPARE_LDAK_SUMMARY_STATISTICS.out.summary_statistics.map { meta, summary_statistics -> [meta.summary_statistics_id, meta, summary_statistics] }

    def ch_sumher_invocations = ch_sumher_requests
        .map { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, tagging_file -> [meta.summary_statistics_id, meta, tagging_file] }
        .combine(ch_prepared_ldak_summaries, by: 0)
        .multiMap { summary_statistics_id, request_meta, tagging_file, _summary_meta, summary_statistics ->
            def route_meta = request_meta + [
                id: summary_statistics_id,
                effective_native_args: getLdakSummaryArguments(request_meta),
            ]
            summary: [route_meta, summary_statistics]
            tagging: [[id: request_meta.reference_bundle_id], tagging_file]
        }

    LDAK_SUMHER(
        ch_sumher_invocations.summary,
        ch_sumher_invocations.tagging,
    )

    def ch_sumcors_left = ch_sumcors_requests
        .map { meta, _hapmap3_snplist, _reference_ld_scores, _regression_weights, tagging_file -> [meta.left_summary_statistics_id, meta, tagging_file] }
        .combine(ch_prepared_ldak_summaries, by: 0)
        .map { _left_summary_statistics_id, request_meta, tagging_file, _left_meta, left_summary_statistics ->
            [request_meta.right_summary_statistics_id, request_meta, tagging_file, left_summary_statistics]
        }

    def ch_sumcors_invocations = ch_sumcors_left
        .combine(ch_prepared_ldak_summaries, by: 0)
        .multiMap { _right_summary_statistics_id, request_meta, tagging_file, left_summary_statistics, right_meta, right_summary_statistics ->
            def route_meta = request_meta + [
                id: request_meta.left_summary_statistics_id,
                effective_native_args: getLdakSummaryArguments(request_meta),
            ]
            left: [route_meta, left_summary_statistics]
            right: [right_meta, right_summary_statistics]
            tagging: [[id: request_meta.reference_bundle_id], tagging_file]
        }

    LDAK_SUMCORS(
        ch_sumcors_invocations.left,
        ch_sumcors_invocations.right,
        ch_sumcors_invocations.tagging,
    )

    emit:
    sumher_hers                    = LDAK_SUMHER.out.hers // channel: [ val(meta), path(native.hers) ], one per ldak_sumher request
    sumher_categories              = LDAK_SUMHER.out.categories // channel: [ val(meta), path(native.cats) ]
    sumher_shares                  = LDAK_SUMHER.out.shares // channel: [ val(meta), path(native.share) ]
    sumher_enrichments             = LDAK_SUMHER.out.enrichments // channel: [ val(meta), path(native.enrich) ]
    sumher_extra                   = LDAK_SUMHER.out.extra // channel: [ val(meta), path(native.extra) ]
    sumher_cross                   = LDAK_SUMHER.out.cross // channel: [ val(meta), path(native.cross) ]
    sumher_taus                    = LDAK_SUMHER.out.taus // channel: [ val(meta), path(native.taus) ]
    sumher_labels                  = LDAK_SUMHER.out.labels // channel: [ val(meta), path(native.labels) ]
    sumher_progress                = LDAK_SUMHER.out.progress // channel: [ val(meta), path(native.progress) ]
    sumher_overlap                 = LDAK_SUMHER.out.overlap // channel: [ val(meta), path(native.overlap) ]
    sumher_hers_liability          = LDAK_SUMHER.out.hers_liability // channel: [ val(meta), path(native.hers.liab) ], optional
    sumher_categories_liability    = LDAK_SUMHER.out.categories_liability // channel: [ val(meta), path(native.cats.liab) ], optional
    sumher_liability_factor        = LDAK_SUMHER.out.liability_factor // channel: [ val(meta), path(native.factor) ], optional
    sumher_log                     = LDAK_SUMHER.out.log // channel: [ val(meta), path(native.log) ]
    sumcors_correlations           = LDAK_SUMCORS.out.correlations // channel: [ val(meta), val(meta2), path(native.cors) ]
    sumcors_correlations_full      = LDAK_SUMCORS.out.correlations_full // channel: [ val(meta), val(meta2), path(native.cors.full) ]
    sumcors_labels                 = LDAK_SUMCORS.out.labels // channel: [ val(meta), val(meta2), path(native.labels) ]
    sumcors_progress               = LDAK_SUMCORS.out.progress // channel: [ val(meta), val(meta2), path(native.progress) ]
    sumcors_overlap                = LDAK_SUMCORS.out.overlap // channel: [ val(meta), val(meta2), path(native.overlap) ]
    sumcors_correlations_liability = LDAK_SUMCORS.out.correlations_liability // channel: [ val(meta), val(meta2), path(native.cors.liab) ], optional
    sumcors_log                    = LDAK_SUMCORS.out.log // channel: [ val(meta), val(meta2), path(native.log) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// The effective native argument list for one summary-scale LDAK request. A request's own `--cutoff` or
// `--truncate` is honoured; otherwise the pipeline's default rare-variant cutoff applies. Liability-scale
// conversion is requested only when the declaring endpoints supply a complete prevalence pair, and SumCors
// emits the declared left-then-right order.
def getLdakSummaryArguments(meta) {
    def args = new ArrayList(meta.native_args)
    def option_names = args.findAll { token -> token instanceof String && token.startsWith('--') }.collect { token -> token.split('=', 2)[0] }
    if (!option_names.contains('--cutoff') && !option_names.contains('--truncate')) {
        args.addAll(['--cutoff', '0.01'])
    }
    if (meta.method == 'ldak_sumher' && meta.is_binary && meta.population_prevalence != null && meta.sample_prevalence != null) {
        args.addAll(['--prevalence', meta.population_prevalence.toString(), '--ascertainment', meta.sample_prevalence.toString()])
    }
    if (meta.method == 'ldak_sumcors' && meta.left_is_binary && meta.right_is_binary && meta.left_population_prevalence != null && meta.left_sample_prevalence != null && meta.right_population_prevalence != null && meta.right_sample_prevalence != null) {
        args.addAll(
            [
                '--prevalence',
                meta.left_population_prevalence.toString(),
                '--ascertainment',
                meta.left_sample_prevalence.toString(),
                '--prevalence2',
                meta.right_population_prevalence.toString(),
                '--ascertainment2',
                meta.right_sample_prevalence.toString(),
            ]
        )
    }
    return args
}
