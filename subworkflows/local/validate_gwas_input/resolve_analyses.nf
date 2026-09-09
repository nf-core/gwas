include { getInternalSummaryMetadata ; getSummaryStatisticsId } from './identity_helpers'
include {
    getAnalysisSettings ;
    getMethodRoutes ;
    normaliseCellValue ;
    tokenizeMethodSelector ;
    validateHeritabilityTraitSupport ;
    validateMethodConditionedColumns ;
    validateMethodSelectors ;
    validateTraitColumns
} from './manifest_contracts'
include { getMethodCapabilities            } from './method_registry'

def getPrevalenceConsumers(relationship_rows, summary_statistics_rows) {
    def pair_prevalence_analysis_ids = relationship_rows.findAll { row ->
        tokenizeMethodSelector(row[0].relationship_methods).any { method ->
            def capability = getMethodCapabilities()[method]
            capability && capability.domain == 'pairwise' && capability.prevalence.population != 'not_consumed'
        }
    }.collectMany { row -> [normaliseCellValue(row[0].left_analysis_id), normaliseCellValue(row[0].right_analysis_id)] }.findAll { analysis_id -> analysis_id }.collect { analysis_id -> analysis_id.toString() } as Set
    def pair_prevalence_summary_statistics_ids = relationship_rows.findAll { row ->
        tokenizeMethodSelector(row[0].relationship_methods).any { method ->
            def capability = getMethodCapabilities()[method]
            capability && capability.domain == 'pairwise' && capability.prevalence.population != 'not_consumed'
        }
    }.collectMany { row -> [normaliseCellValue(row[0].left_summary_statistics_id), normaliseCellValue(row[0].right_summary_statistics_id)] }.findAll { summary_statistics_id -> summary_statistics_id }.collect { summary_statistics_id -> summary_statistics_id.toString() } as Set
    def summary_prevalence_analysis_ids = summary_statistics_rows.findAll { row ->
        def summary_statistics_id = normaliseCellValue(row[0].id)?.toString()
        normaliseCellValue(row[0].producer_analysis_id) && (tokenizeMethodSelector(row[0].heritability_methods).any { method ->
            def capability = getMethodCapabilities()[method]
            capability && capability.prevalence.population != 'not_consumed'
        } || summary_statistics_id in pair_prevalence_summary_statistics_ids)
    }.collect { row -> normaliseCellValue(row[0].producer_analysis_id).toString() } as Set
    return [
        pair_analysis_ids: pair_prevalence_analysis_ids,
        summary_analysis_ids: summary_prevalence_analysis_ids,
    ]
}

def resolveAnalyses(analysis_rows, analysis_columns, analysis_manifest, cohort_manifest, cohorts_by_id, relationship_rows, prevalence_consumers, method_options_by_analysis, errors) {
    def analyses_by_id = [:]
    def line_by_analysis_id = [:]
    def validated_analyses = []
    def generated_summaries_by_id = [:]
    def pair_prevalence_analysis_ids = prevalence_consumers.pair_analysis_ids
    def summary_prevalence_analysis_ids = prevalence_consumers.summary_analysis_ids
    analysis_rows.eachWithIndex { row, index ->
        def line = index + 2
        def analysis_meta = row[0]
        def cells = [analysis_columns, row[1..-1]].transpose().collectEntries()
        def analysis_id = analysis_meta.id
        def cohort_id = analysis_meta.cohort
        def reject = { field, message ->
            def named = field instanceof List ? field : [field]
            def label = named.size() > 1
                ? "fields ${named.collect { name -> "'${name}'" }.join(', ')}"
                : "field '${named.first()}'"
            errors << "  - ${analysis_manifest} row ${line} (analysis_id '${analysis_id}'), ${label}: ${message}"
        }

        if (line_by_analysis_id.containsKey(analysis_id)) {
            reject.call('analysis_id', "duplicate analysis_id '${analysis_id}', already declared on row ${line_by_analysis_id[analysis_id]}")
        }
        else {
            line_by_analysis_id[analysis_id] = line
        }
        def cohort = cohorts_by_id[cohort_id]
        if (!cohort) {
            reject.call('cohort_id', "undefined cohort_id '${cohort_id}', not declared in cohort manifest '${cohort_manifest}'")
        }

        def association_methods = tokenizeMethodSelector(analysis_meta.association_methods)
        def heritability_methods = tokenizeMethodSelector(analysis_meta.heritability_methods)
        validateMethodSelectors(association_methods, heritability_methods, reject, relationship_rows ? true : false)
        def routes = getMethodRoutes(association_methods, heritability_methods)
        def settings = getAnalysisSettings(analysis_meta)
        def is_binary = analysis_meta.trait_type == 'binary'
        validateHeritabilityTraitSupport(is_binary, heritability_methods, reject)
        validateTraitColumns(is_binary, settings, reject)
        validateMethodConditionedColumns(
            settings,
            routes,
            reject,
            analysis_id in pair_prevalence_analysis_ids || analysis_id in summary_prevalence_analysis_ids,
        )

        if (cohort) {
            // REGENIE needs the user's categorical column names to expand the combined covariate file.
            def cat_covariate_names = []
            if (cells.cat_covariates) {
                def header = file(cells.cat_covariates).readLines().find { header_line -> header_line.trim() }
                cat_covariate_names = (header.contains('\t') ? header.split('\t', -1) : header.trim().split(/\s+/)).drop(2)
            }
            // The declared genotype view identity is added as a key only when the cohort declared one, and is
            // never added as a `null`. Analysis metadata is a task-hash input everywhere downstream, so an
            // unconditional `genotype_view_id: null` would change every existing analysis's hash and every
            // pinned snapshot for cohorts that declare nothing.
            def resolved_meta = analysis_meta + [
                build: cohort.meta.build,
                ancestry: cohort.meta.ancestry,
                association_methods: association_methods,
                heritability_methods: heritability_methods,
                genotype_format: cohort.genotype_format,
                is_binary: is_binary,
                cat_covariate_names: cat_covariate_names,
                case_value: settings.case_value == null ? null : settings.case_value.toString(),
                control_value: settings.control_value == null ? null : settings.control_value.toString(),
                population_prevalence: settings.population_prevalence,
                sample_prevalence: settings.sample_prevalence,
                method_options: method_options_by_analysis[analysis_id],
            ] + (cohort.genotype_view_id ? [genotype_view_id: cohort.genotype_view_id] : [:])
            def validated = [
                resolved_meta,
                cohort.genotype_files,
                cells.phenotype,
                cells.quant_covariates,
                cells.cat_covariates,
                method_options_by_analysis[analysis_id].ldak.predictor_extract,
                method_options_by_analysis[analysis_id].ldak.weights,
            ]
            validated_analyses << validated
            analyses_by_id[analysis_id] = validated
            association_methods.each { association_method ->
                def summary_statistics_id = getSummaryStatisticsId(analysis_id, association_method)
                generated_summaries_by_id[summary_statistics_id] = getInternalSummaryMetadata(resolved_meta, association_method)
            }
        }
    }
    return [
        analyses: validated_analyses,
        analyses_by_id: analyses_by_id,
        generated_summaries_by_id: generated_summaries_by_id,
        line_by_analysis_id: line_by_analysis_id,
    ]
}
