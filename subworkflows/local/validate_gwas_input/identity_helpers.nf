def getSummaryStatisticsId(analysis_id, association_method) {
    return "${analysis_id}--${association_method}".toString()
}

def getInternalSummaryMetadata(meta, association_method) {
    def summary_statistics_id = getSummaryStatisticsId(meta.id, association_method)
    return [
        id: summary_statistics_id,
        summary_statistics_id: summary_statistics_id,
        trait: meta.trait,
        trait_id: meta.trait,
        trait_type: meta.trait_type,
        is_binary: meta.is_binary,
        population_prevalence: meta.population_prevalence,
        sample_prevalence: meta.sample_prevalence,
        build: meta.build,
        ancestry: meta.ancestry,
        source_kind: 'pipeline_generated',
        source_format: "pipeline_${association_method}",
        source_method: association_method,
        source_release: null,
        source_name: "${meta.id}.${association_method}",
        producer_analysis_id: meta.id,
        producer_association_method: association_method,
        access_constraints: null,
    ]
}
