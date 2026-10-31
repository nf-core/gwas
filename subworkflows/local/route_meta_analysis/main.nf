// Pipeline-owned ordered membership, model selection, result assembly and derived-summary identity.

// MODULE: Local to the pipeline
include { GWASLAB_METAANALYZE          } from '../../../modules/local/gwaslab/metaanalyze/main'
include { METASOFT_RE2                 } from '../../../modules/local/metasoft/re2/main'
include { FINALISE_META_ANALYSIS       } from '../../../modules/local/finalise_meta_analysis/main'

// SUBWORKFLOW: Reusable composition of local components
include { SUMSTATS_METAANALYSIS_MRMEGA } from '../sumstats_metaanalysis_mrmega/main'

workflow ROUTE_META_ANALYSIS {
    take:
    ch_requests // channel: [ val(meta), val(ordered_source_summary_statistics_ids) ]
    ch_summary_statistics // channel: [ val(meta), path(gwaslab_summary_statistics) ], base summaries only

    main:
    def ch_parents = ch_requests
        .flatMap { meta, source_ids ->
            source_ids
                .withIndex()
                .collect { source_id, ordinal ->
                    [source_id, groupKey(meta, source_ids.size()), ordinal]
                }
        }
        .combine(
            ch_summary_statistics.map { meta, source -> [meta.summary_statistics_id, source] },
            by: 0
        )
        .map { source_id, request_key, ordinal, source -> [request_key, ordinal, source_id, source] }
        .groupTuple(remainder: true)
        .map { request_key, ordinals, source_ids, sources ->
            def meta = request_key.getGroupTarget()
            def ordered = orderMetaAnalysisParents(meta, ordinals, source_ids, sources)
            [meta, ordered.sources, ordered.ids]
        }

    // Every request yields the fixed-effects table and the aligned adapter inputs the optional programmes read.
    // `random` is GWASLab's own optional model and is selected on the call; `re2` and `mrmega` select the two
    // native programmes below. Model selection is pipeline policy and stays in this route.
    GWASLAB_METAANALYZE(
        ch_parents.map { meta, sources, source_ids ->
            [meta, sources, source_ids, 'gwaslab', meta.build, meta.meta_analysis_models.contains('random')]
        }
    )

    METASOFT_RE2(
        GWASLAB_METAANALYZE.out.metasoft_input.filter { meta, _matrix -> meta.meta_analysis_models.contains('re2') }
    )

    // The ordered parent identities are the study names of the aligned columns: the gather above fails the run
    // unless they equal the request's declared `source_summary_statistics_ids`.
    SUMSTATS_METAANALYSIS_MRMEGA(
        GWASLAB_METAANALYZE.out.mrmega_input
            .filter { meta, _aligned -> meta.meta_analysis_models.contains('mrmega') }
            .map { meta, aligned -> [meta, aligned, meta.source_summary_statistics_ids, meta.trait_type] }
    )

    // Each unselected model supplies one empty slot per request. A selected model supplies its native file;
    // keyed joins enforce this route's one-result-per-selected-model contract without inspecting statistics.
    def ch_random = ch_requests
        .filter { meta, _source_ids -> !meta.meta_analysis_models.contains('random') }
        .map { meta, _source_ids -> [meta.request_id, []] }
        .mix(GWASLAB_METAANALYZE.out.random_effects.map { meta, result -> [meta.request_id, result] })
    def ch_re2 = ch_requests
        .filter { meta, _source_ids -> !meta.meta_analysis_models.contains('re2') }
        .map { meta, _source_ids -> [meta.request_id, []] }
        .mix(METASOFT_RE2.out.result.map { meta, result -> [meta.request_id, result] })
    def ch_mrmega = ch_requests
        .filter { meta, _source_ids -> !meta.meta_analysis_models.contains('mrmega') }
        .map { meta, _source_ids -> [meta.request_id, []] }
        .mix(SUMSTATS_METAANALYSIS_MRMEGA.out.result.map { meta, result -> [meta.request_id, result] })

    FINALISE_META_ANALYSIS(
        GWASLAB_METAANALYZE.out.fixed.map { meta, fixed -> [meta.request_id, meta, fixed] }.join(ch_random, failOnDuplicate: true, failOnMismatch: true).join(ch_re2, failOnDuplicate: true, failOnMismatch: true).join(ch_mrmega, failOnDuplicate: true, failOnMismatch: true).map { _request_id, meta, fixed, random, re2, mrmega -> [meta, fixed, random, re2, mrmega] }
    )

    emit:
    summary_statistics = FINALISE_META_ANALYSIS.out.summary_statistics // channel: [ val(meta), path(derived_gwaslab_summary) ]
    fixed              = GWASLAB_METAANALYZE.out.fixed // channel: [ val(meta), path(native_fixed) ]
    random_effects     = GWASLAB_METAANALYZE.out.random_effects // channel: [ val(meta), path(native_random) ]
    gwaslab_log        = GWASLAB_METAANALYZE.out.log // channel: [ val(meta), path(native_log) ]
    metasoft_result    = METASOFT_RE2.out.result // channel: [ val(meta), path(native_result) ]
    metasoft_log       = METASOFT_RE2.out.log // channel: [ val(meta), path(native_log) ]
    mrmega_result      = SUMSTATS_METAANALYSIS_MRMEGA.out.result // channel: [ val(meta), path(native_result) ]
    mrmega_log         = SUMSTATS_METAANALYSIS_MRMEGA.out.mrmega_log // channel: [ val(meta), path(native_log) ]
}

def orderMetaAnalysisParents(meta, ordinals, source_ids, sources) {
    def expected = meta.source_summary_statistics_ids
    def ordered = [ordinals, source_ids, sources].transpose().sort { left, right -> left[0] <=> right[0] }
    if (ordered.collect { row -> row[1] } != expected) {
        error("[nf-core/gwas] ERROR: Meta-analysis request '${meta.request_id}' did not resolve exactly its declared ordered summary parents")
    }
    return [ids: ordered.collect { row -> row[1] }, sources: ordered.collect { row -> row[2] }]
}
