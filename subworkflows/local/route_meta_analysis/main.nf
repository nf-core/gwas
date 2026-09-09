// Pipeline-owned ordered membership, result assembly and derived-summary identity.
include { COMMON_VARIANT_META_ANALYSIS } from '../common_variant_meta_analysis/main'
include { FINALISE_META_ANALYSIS } from '../../../modules/local/finalise_meta_analysis/main'

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
            [meta, ordered.sources, ordered.ids, 'gwaslab', meta.build, meta.trait_type, meta.meta_analysis_models, meta.axes]
        }

    COMMON_VARIANT_META_ANALYSIS(ch_parents)

    // Each unselected model supplies one empty slot per request. A selected model supplies its native file;
    // keyed joins enforce this route's one-result-per-selected-model contract without inspecting statistics.
    def ch_random = ch_requests
        .filter { meta, _source_ids -> !meta.meta_analysis_models.contains('random') }
        .map { meta, _source_ids -> [meta.request_id, []] }
        .mix(COMMON_VARIANT_META_ANALYSIS.out.random_effects.map { meta, result -> [meta.request_id, result] })
    def ch_re2 = ch_requests
        .filter { meta, _source_ids -> !meta.meta_analysis_models.contains('re2') }
        .map { meta, _source_ids -> [meta.request_id, []] }
        .mix(COMMON_VARIANT_META_ANALYSIS.out.metasoft_result.map { meta, result -> [meta.request_id, result] })
    def ch_mrmega = ch_requests
        .filter { meta, _source_ids -> !meta.meta_analysis_models.contains('mrmega') }
        .map { meta, _source_ids -> [meta.request_id, []] }
        .mix(COMMON_VARIANT_META_ANALYSIS.out.mrmega_result.map { meta, result -> [meta.request_id, result] })

    FINALISE_META_ANALYSIS(
        COMMON_VARIANT_META_ANALYSIS.out.fixed.map { meta, fixed -> [meta.request_id, meta, fixed] }.join(ch_random, failOnDuplicate: true, failOnMismatch: true).join(ch_re2, failOnDuplicate: true, failOnMismatch: true).join(ch_mrmega, failOnDuplicate: true, failOnMismatch: true).map { _request_id, meta, fixed, random, re2, mrmega -> [meta, fixed, random, re2, mrmega] }
    )

    emit:
    summary_statistics = FINALISE_META_ANALYSIS.out.summary_statistics // channel: [ val(meta), path(derived_gwaslab_summary) ]
    fixed = COMMON_VARIANT_META_ANALYSIS.out.fixed // channel: [ val(meta), path(native_fixed) ]
    random_effects = COMMON_VARIANT_META_ANALYSIS.out.random_effects // channel: [ val(meta), path(native_random) ]
    gwaslab_log = COMMON_VARIANT_META_ANALYSIS.out.gwaslab_log // channel: [ val(meta), path(native_log) ]
    metasoft_result = COMMON_VARIANT_META_ANALYSIS.out.metasoft_result // channel: [ val(meta), path(native_result) ]
    metasoft_log = COMMON_VARIANT_META_ANALYSIS.out.metasoft_log // channel: [ val(meta), path(native_log) ]
    mrmega_result = COMMON_VARIANT_META_ANALYSIS.out.mrmega_result // channel: [ val(meta), path(native_result) ]
    mrmega_log = COMMON_VARIANT_META_ANALYSIS.out.mrmega_log // channel: [ val(meta), path(native_log) ]
}

def orderMetaAnalysisParents(meta, ordinals, source_ids, sources) {
    def expected = meta.source_summary_statistics_ids
    def ordered = [ordinals, source_ids, sources].transpose().sort { left, right -> left[0] <=> right[0] }
    if (ordered.collect { row -> row[1] } != expected) {
        error("[nf-core/gwas] ERROR: Meta-analysis request '${meta.request_id}' did not resolve exactly its declared ordered summary parents")
    }
    return [ids: ordered.collect { row -> row[1] }, sources: ordered.collect { row -> row[2] }]
}
