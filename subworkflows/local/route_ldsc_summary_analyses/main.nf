// Route nf-core/gwas standalone CBIIT Python 3 LDSC requests — heritability (H2) and genetic correlation
// (RG) — from GWASLab-standard summary statistics through one shared content-addressed munging to native
// LDSC logs. H2 and RG are one controller because they share that munging: a GWASLab summary
// consumed by a unary request and by either side of any number of pair requests is munged exactly once.
// This is pipeline routing, reuse identity and scientific/publication policy, not an nf-core/modules
// submission candidate. Every constituent process reports directly to the run-wide versions topic, so this
// subworkflow emits no versions, and it reads no params, no workflow and no projectDir.

// MODULES: Upstream-ready components used inside a pipeline-local route
include { LDSC_H2            } from '../../../modules/local/ldsc/h2/main'
include { LDSC_MUNGESUMSTATS } from '../../../modules/local/ldsc/mungesumstats/main'
include { LDSC_RG            } from '../../../modules/local/ldsc/rg/main'

// FUNCTION: Local to the pipeline
include { digestFileBytes    } from '../utils_nfcore_gwas_pipeline'
include { digestIdentityText } from '../utils_nfcore_gwas_pipeline'

workflow ROUTE_LDSC_SUMMARY_ANALYSES {
    take:
    ch_h2_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ], the ldsc_h2 unary requests; an LDSC-family bundle supplies [] for the tagging file
    ch_rg_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ], the ldsc_rg pair requests, same convention
    ch_summary_statistics // channel: [ val(meta), path(gwaslab_summary_statistics) ], one element per distinct summary_statistics_id

    main:

    // Munging belongs to a GWASLab summary plus the exact HapMap3 allele-universe bytes, not to a
    // downstream request or its regression reference/weights. A content-derived key therefore lets unary,
    // pairwise, primary and named sensitivity requests reuse the same expensive preparation without making
    // ancestry, bundle names, H2/RG native arguments or output identity part of that derivation.
    def ch_ldsc_munging_requests = ch_h2_requests
        .map { meta, hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file ->
            def key = getLdscMungingKey(meta.summary_statistics_id, hapmap3_snplist)
            [meta.summary_statistics_id, key, hapmap3_snplist]
        }
        .mix(
            ch_rg_requests.flatMap { meta, hapmap3_snplist, _reference_ld_scores, _regression_weights, _tagging_file ->
                [meta.left_summary_statistics_id, meta.right_summary_statistics_id].collect { summary_statistics_id ->
                    [summary_statistics_id, getLdscMungingKey(summary_statistics_id, hapmap3_snplist), hapmap3_snplist]
                }
            }
        )
        .unique { _summary_statistics_id, key, _hapmap3_snplist -> key }

    def ch_summary_by_id = ch_summary_statistics.map { meta, summary_statistics -> [meta.summary_statistics_id, meta, summary_statistics] }

    def ch_ldsc_munging_invocations = ch_ldsc_munging_requests
        .combine(ch_summary_by_id, by: 0)
        .multiMap { summary_statistics_id, key, hapmap3_snplist, summary_meta, summary_statistics ->
            def munging_meta = summary_meta + [
                id: key,
                munging_key: key,
                summary_statistics_id: summary_statistics_id,
                hapmap3_sha256: digestFileBytes(hapmap3_snplist),
                munging_adapter_contract: 'gwaslab_to_ldsc_sumstats_v1',
            ]
            sumstats: [munging_meta, summary_statistics]
            merge_alleles: [[id: key], hapmap3_snplist]
        }

    LDSC_MUNGESUMSTATS(
        ch_ldsc_munging_invocations.sumstats,
        ch_ldsc_munging_invocations.merge_alleles,
    )

    def ch_ldsc_munged = LDSC_MUNGESUMSTATS.out.munged_sumstats.map { meta, munged_sumstats -> [meta.munging_key, meta, munged_sumstats] }

    // Each request runs once, with native liability conversion only when the binary endpoint declares
    // both prevalences. Request identity and LD/weight resources are joined after shared munging.
    def ch_ldsc_h2_invocations = ch_h2_requests
        .map { meta, hapmap3_snplist, reference_ld_scores, regression_weights, _tagging_file ->
            [getLdscMungingKey(meta.summary_statistics_id, hapmap3_snplist), meta, reference_ld_scores, regression_weights]
        }
        .combine(ch_ldsc_munged, by: 0)
        .multiMap { key, meta, reference_ld_scores, regression_weights, _munging_meta, munged_sumstats ->
            def liability = meta.is_binary && meta.population_prevalence != null && meta.sample_prevalence != null
            def route_meta = meta + [
                munging_keys: [key],
                native_scale: liability ? 'liability' : 'observed',
                effective_population_prevalence: liability ? [meta.population_prevalence] : [],
                effective_sample_prevalence: liability ? [meta.sample_prevalence] : [],
            ]
            sumstats: [route_meta, munged_sumstats]
            reference_ld_scores: [[id: meta.reference_bundle_id], reference_ld_scores]
            regression_weights: [[id: meta.reference_bundle_id], regression_weights]
        }

    LDSC_H2(
        ch_ldsc_h2_invocations.sumstats,
        ch_ldsc_h2_invocations.reference_ld_scores,
        ch_ldsc_h2_invocations.regression_weights,
    )

    // Pair requests preserve declared left/right order. Both endpoint munging keys are resolved against the
    // one HapMap3 resource selected by this request, then the request-owned LD-score and regression-weight
    // directories are passed unchanged to RG.
    def ch_ldsc_rg_left = ch_rg_requests
        .map { meta, hapmap3_snplist, reference_ld_scores, regression_weights, _tagging_file ->
            def left_key = getLdscMungingKey(meta.left_summary_statistics_id, hapmap3_snplist)
            [left_key, meta, hapmap3_snplist, reference_ld_scores, regression_weights]
        }
        .combine(ch_ldsc_munged, by: 0)
        .map { left_key, meta, hapmap3_snplist, reference_ld_scores, regression_weights, _left_munging_meta, left_sumstats ->
            def right_key = getLdscMungingKey(meta.right_summary_statistics_id, hapmap3_snplist)
            [right_key, left_key, meta, reference_ld_scores, regression_weights, left_sumstats]
        }

    def ch_ldsc_rg_invocations = ch_ldsc_rg_left
        .combine(ch_ldsc_munged, by: 0)
        .multiMap { right_key, left_key, meta, reference_ld_scores, regression_weights, left_sumstats, _right_munging_meta, right_sumstats ->
            def has_binary = meta.left_is_binary || meta.right_is_binary
            def complete = [
                [binary: meta.left_is_binary, population: meta.left_population_prevalence, sample: meta.left_sample_prevalence],
                [binary: meta.right_is_binary, population: meta.right_population_prevalence, sample: meta.right_sample_prevalence],
            ].every { endpoint -> !endpoint.binary || (endpoint.population != null && endpoint.sample != null) }
            def liability = has_binary && complete
            def population = [
                meta.left_is_binary ? meta.left_population_prevalence : 'nan',
                meta.right_is_binary ? meta.right_population_prevalence : 'nan',
            ]
            def sample = [
                meta.left_is_binary ? meta.left_sample_prevalence : 'nan',
                meta.right_is_binary ? meta.right_sample_prevalence : 'nan',
            ]
            def route_meta = meta + [
                munging_keys: [left_key, right_key],
                native_scale: liability ? 'liability' : 'observed',
                effective_population_prevalence: liability ? population : [],
                effective_sample_prevalence: liability ? sample : [],
            ]
            sumstats: [route_meta, left_sumstats, right_sumstats]
            reference_ld_scores: [[id: meta.reference_bundle_id], reference_ld_scores]
            regression_weights: [[id: meta.reference_bundle_id], regression_weights]
        }

    LDSC_RG(
        ch_ldsc_rg_invocations.sumstats,
        ch_ldsc_rg_invocations.reference_ld_scores,
        ch_ldsc_rg_invocations.regression_weights,
    )

    emit:
    h2_log = LDSC_H2.out.log.map { meta, log -> [stripRoutingIdentity(meta), log] } // channel: [ val(meta), path(native.log) ]
    rg_log = LDSC_RG.out.log.map { meta, log -> [stripRoutingIdentity(meta), log] } // channel: [ val(meta), path(native.log) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// The reusable munging identity of one GWASLab summary. It is exactly the adapter contract, the summary
// statistics identity and the SHA-256 of the HapMap3 SNP-list bytes. Request identity, native arguments, LD
// scores, regression weights, reference-bundle names, ancestry labels and sensitivity names are deliberately
// unreachable from here so that they cannot make two scientifically identical mungings distinct.
def getLdscMungingKey(summary_statistics_id, hapmap3_snplist) {
    return digestIdentityText(
        [
            'adapter=gwaslab_to_ldsc_sumstats_v1',
            "summary_statistics_id=${summary_statistics_id}",
            "hapmap3_sha256=${digestFileBytes(hapmap3_snplist)}",
        ].join('\n')
    )
}

// Drop the controller's private reuse identity from a result record. The focal scientific identity — request,
// method, relationship and endpoint attribution, trait identity and prevalence declarations — is retained.
def stripRoutingIdentity(meta) {
    return meta.findAll { key, _value -> key != 'munging_keys' }
}
