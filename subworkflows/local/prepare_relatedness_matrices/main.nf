// Build each relatedness base artifact once, derive each requested child once, and restore focal metadata.
// Every component reports on the run-wide versions topic, so this subworkflow emits no versions.

// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
include { PLINK_PREPARE_GRM_GCTA           } from '../plink_prepare_grm_gcta/main'
include { PLINK_PREPARE_LDMS_PLAN_GCTA     } from '../plink_prepare_ldms_plan_gcta/main'
include { PLINK_PREPARE_GRM_LDMS_GCTA      } from '../plink_prepare_grm_ldms_gcta/main'
include { PLINK_PREPARE_GRM_LDAK           } from '../plink_prepare_grm_ldak/main'
include { PLINK_PREPARE_GRM_MPH            } from '../plink_prepare_grm_mph/main'
include { PLINK_PREPARE_GRM_LDMS_MPH       } from '../plink_prepare_grm_ldms_mph/main'

// MODULE: Local to the pipeline
include { GCTA_MAKEBKSPARSE                } from '../../../modules/local/gcta/makebksparse/main'

// FUNCTION: Local to the pipeline
include { buildScientificArtifactKey       } from '../utils_nfcore_gwas_pipeline'
include { canonicaliseScientificIdentifier } from '../utils_nfcore_gwas_pipeline'
include { canonicaliseScientificValue      } from '../utils_nfcore_gwas_pipeline'
include { digestFileBytes                  } from '../utils_nfcore_gwas_pipeline'
include { getMatrixKindContract            } from '../validate_gwas_input/method_registry'
include { getMethodCapabilities            } from '../validate_gwas_input/method_registry'

workflow PREPARE_RELATEDNESS_MATRICES {
    take:
    ch_analyses // channel: [ val(meta), [ path(genotype_file), ... ], path(ldak_weights) ], [] when absent
    ch_cohort_native_genotypes // channel: [ val(cohort_meta), path(primary_genotype), path(variant_file), path(sample_file) ], format-polymorphic
    ch_plink1_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam) ], analyses needing PLINK 1
    ch_cohort_native_view_keys // channel: [ val(cohort_id), val(native_view_key) ], one per cohort
    ch_cohort_plink1_view_keys // channel: [ val(cohort_id), val(plink1_view_key) ], only cohorts with a PLINK 1 view
    gcta_grm_parts // channel: val(gcta_grm_parts), run/profile GCTA GRM partition count

    main:

    // The highest integer chromosome code the MPH matrix families count as autosomal, matching GCTA's
    // `--autosome-num` default. It is a constant here rather than a pipeline parameter because it is not a
    // per-analysis choice: GCTA restricts `--make-grm` to the autosomes and refuses `--ld-score-region`
    // outright on anything else, while MPH applies no chromosome rule at all and reports nothing, so the two
    // engines would otherwise estimate over different predictor sets from one bundle. The universe travels
    // into the `mph_dense` reuse key symbolically rather than numerically, so changing this constant cannot
    // silently reuse today's bytes under today's key.
    def autosome_count = 22

    // Final consumer requests preserve focal metadata beside an explicit base/derived artifact graph. The
    // view-compatibility key a request is built against is chosen here, at the caller, from the representation
    // the kind's builder actually reads: a `plink1` kind is keyed by the cohort's PLINK 1 view, everything
    // else by its native view. Selecting per kind is what lets the key channel that is total for that kind be
    // the only one it is combined with — a `plink1` kind on a cohort with no PLINK 1 view cannot arise,
    // because the same registry predicate produces the demand and the kind.
    def ch_kind_requests = ch_analyses.flatMap { meta, _genotype_files, ldak_weights ->
        getRelatednessMatrixKinds(meta).collect { kind -> [meta.cohort, meta, kind, ldak_weights] }
    }

    def ch_kind_bundles = ch_kind_requests.branch { _cohort_id, _meta, kind, _ldak_weights ->
        plink1: getMatrixKindContract()[kind].genotype_bundle == 'plink1'
        native: true
    }

    def ch_requests = ch_kind_bundles.plink1
        .combine(ch_cohort_plink1_view_keys, by: 0)
        .mix(ch_kind_bundles.native.combine(ch_cohort_native_view_keys, by: 0))
        .map { _cohort_id, meta, kind, ldak_weights, view_compatibility_key ->
            def weights_policy = kind == 'ldak_kinship' ? meta.method_options.ldak.weights_policy : 'equal'
            def weights_identity = kind == 'ldak_kinship' ? getLdakWeightsIdentity(ldak_weights, weights_policy) : [mode: 'equal']
            def gcta_extract = kind in ['gcta_dense', 'gcta_sparse'] && !meta.relationship_id ? meta.method_options.gcta.grm_extract : []
            def extract_identity = kind in ['gcta_dense', 'gcta_sparse'] ? getMethodResourceIdentity(gcta_extract) : [mode: 'all']
            def request = buildRelatednessArtifactRequest(meta, view_compatibility_key, kind, weights_identity, extract_identity, gcta_extract)
            def weights_file = kind == 'ldak_kinship' ? ldak_weights ?: [] : []
            [request.selected_key, meta, request, weights_file]
        }

    // Base construction is reduced independently of final requested form. Sparse cutoffs and LDAK filtering
    // cannot fragment these records because they occur only on child requests.
    def ch_base_requests = ch_requests
        .map { _selected_key, _meta, request, weights_file -> [request.base.key, request.base, weights_file] }
        .unique { key, _base, _weights_file -> key }
        .map { key, base, weights_file ->
            def matrix_meta = buildRelatednessArtifactMeta(base)
            [key, matrix_meta, base, weights_file, gcta_grm_parts]
        }

    // Genotypes are paired to a base request by the view identity the request was keyed against, not by
    // cohort. Two `cohort_id`s over byte-identical files now collapse to one base key at the `unique` above,
    // so a cohort-keyed pairing would attribute the surviving shared artifact — its tag, its staged filenames
    // and its trace row — to whichever request the flatMap happened to emit first. The published matrix
    // directory is already named by the key, so the identity is what the pairing should use.
    //
    // Collapsing several cohorts onto one view therefore has to choose one of their bundles, and the choice
    // must not be first-to-arrive. The bundles are byte-identical but they are distinct files — a PGEN or VCF
    // cohort is projected and imported under its own cohort id, so two cohorts sharing a view produce two
    // parallel tasks writing differently named outputs — and Nextflow hashes a path input by its source path.
    // A first-arrival `unique` would therefore give the shared matrix a different task hash depending on
    // which producer finished first, so an unchanged `-resume` could rebuild every shared matrix, and the
    // GCTA manifest stem would name an arbitrary cohort. Grouping and taking the lowest cohort id makes the
    // survivor a function of the manifest instead of of the scheduler. The cost is that the reduction waits
    // for every cohort's preparation to complete rather than streaming; preparation is one short task per
    // cohort and they run in parallel, so this delays the first matrix build by at most the slowest single
    // preparation, and never serialises them.
    def ch_native_by_view_key = ch_cohort_native_genotypes
        .map { cohort_meta, primary, variant_file, sample_file -> [cohort_meta.id, primary, variant_file, sample_file] }
        .combine(ch_cohort_native_view_keys, by: 0)
        .map { cohort_id, primary, variant_file, sample_file, view_key -> [view_key, cohort_id, primary, variant_file, sample_file] }
        .groupTuple(by: 0)
        .map { view_key, cohort_ids, primaries, variant_files, sample_files ->
            [view_key] + selectLowestCohort(cohort_ids, [primaries, variant_files, sample_files])
        }

    // One PLINK 1 view exists per view identity even when several cohorts or consumers requested it. Several
    // requests on one cohort reach this stream as identical elements, so they are reduced per cohort first
    // and only distinct cohorts are ever compared.
    def ch_plink1_by_view_key = ch_plink1_genotypes
        .map { meta, bed, bim, fam -> [meta.cohort, bed, bim, fam] }
        .unique { cohort_id, _bed, _bim, _fam -> cohort_id }
        .combine(ch_cohort_plink1_view_keys, by: 0)
        .map { cohort_id, bed, bim, fam, view_key -> [view_key, cohort_id, bed, bim, fam] }
        .groupTuple(by: 0)
        .map { view_key, cohort_ids, beds, bims, fams -> [view_key] + selectLowestCohort(cohort_ids, [beds, bims, fams]) }

    // Branch by base type before pairing. `gcta_ldms` reads PLINK 1 and `gcta_dense` reads the native bundle,
    // so the two must not be routed through one combined stream: doing that used to force the LDMS builds to
    // take a native bundle they immediately discarded. The five branches enumerate every base type
    // `getMatrixKindContract()` can produce — `gcta_sparse` has no branch because it is a derivative of the
    // `gcta_dense` base rather than a base of its own.
    def ch_base_by_type = ch_base_requests.branch { _key, _matrix_meta, base, _weights_file, _parts ->
        gcta_dense: base.type == 'gcta_dense'
        gcta_ldms: base.type == 'gcta_ldms'
        ldak_kinship: base.type == 'ldak_kinship'
        mph_dense: base.type == 'mph_dense'
        mph_ldms: base.type == 'mph_ldms'
    }

    def ch_dense_genotypes = ch_base_by_type.gcta_dense
        .map { _key, matrix_meta, base, _weights_file, parts -> [base.view_compatibility_key, matrix_meta, base, parts] }
        .combine(ch_native_by_view_key, by: 0)
        .map { _view_key, matrix_meta, base, parts, primary, variant_file, sample_file -> [matrix_meta, base, parts, primary, variant_file, sample_file] }

    // GCTA resolves one manifest entry from the staged bundle prefix, and `GCTA_MAKEGRMPART` selects
    // `--mbfile` or `--mpfile` from that primary member's extension, so one manifest shape serves both
    // representations. Matching companion basenames remain a native caller contract, enforced at ingress; the
    // staged basename is not part of the scientific base key.
    def ch_dense_manifests = ch_dense_genotypes
        .map { matrix_meta, _base, _parts, primary, _variant_file, _sample_file -> [matrix_meta, primary.baseName] }
        .collectFile { matrix_meta, stem -> ["${matrix_meta.id}.mfile", "${stem}\n"] }
        .map { manifest -> [manifest.baseName, manifest] }

    def ch_dense_inputs = ch_dense_genotypes
        .map { matrix_meta, base, parts, primary, variant_file, sample_file -> [matrix_meta.id, matrix_meta, base, parts, primary, variant_file, sample_file] }
        .join(ch_dense_manifests, failOnDuplicate: true, failOnMismatch: true)
        .multiMap { _matrix_id, matrix_meta, base, parts, primary, variant_file, sample_file, manifest ->
            genotypes: [matrix_meta, manifest, primary, variant_file, sample_file]
            snp_group: [matrix_meta, base.gcta_extract ?: []]
            n_parts: [matrix_meta, parts]
        }

    PLINK_PREPARE_GRM_GCTA(
        ch_dense_inputs.genotypes,
        ch_dense_inputs.snp_group,
        ch_dense_inputs.n_parts,
    )

    // Each sparse child is keyed by its dense parent plus cutoff. Several cutoffs combine with the same
    // completed parent without creating another dense construction.
    def ch_sparse_derivatives = ch_requests
        .filter { _selected_key, _meta, request, _weights_file -> request.derived && request.derived.type == 'gcta_sparse' }
        .map { _selected_key, _meta, request, _weights_file -> [request.derived.key, request.derived] }
        .unique { key, _derived -> key }

    def ch_sparse_inputs = ch_sparse_derivatives
        .map { _key, derived -> [derived.parent_key, derived] }
        .combine(
            PLINK_PREPARE_GRM_GCTA.out.grm_files.map { matrix_meta, grm_files -> [matrix_meta.key, grm_files] },
            by: 0
        )
        .multiMap { _parent_key, derived, grm_files ->
            def matrix_meta = buildRelatednessArtifactMeta(derived)
            grm: [matrix_meta, grm_files]
            cutoff: derived.settings.cutoff
        }

    GCTA_MAKEBKSPARSE(
        ch_sparse_inputs.grm,
        ch_sparse_inputs.cutoff,
    )

    // The LD-by-MAF component plan is its own reusable artifact, one layer above the matrix families that
    // consume it. Two engines now build stratified matrices from one plan, so an identical plan is computed
    // once and its ordered manifest and SNP group files are the only thing either builder reads. The plan is
    // keyed by the same cohort view as the matrix bases plus the three plan settings, so a request that
    // changes any of them gets its own plan and everything downstream of it.
    def ch_plan_requests = ch_base_requests
        .filter { _key, _matrix_meta, base, _weights_file, _parts -> base.plan }
        .map { _key, _matrix_meta, base, _weights_file, _parts -> [base.plan.key, base.plan] }
        .unique { key, _plan -> key }
        .map { _key, plan -> [plan.view_compatibility_key, buildLdmsPlanMeta(plan)] }
        .combine(ch_plink1_by_view_key, by: 0)
        .multiMap { _view_key, plan_meta, bed, bim, fam ->
            genotypes: [plan_meta, bed, bim, fam]
            region_kb: [plan_meta, plan_meta.settings.ld_score_region_kb]
            ld_bins: [plan_meta, plan_meta.settings.ld_bins]
            maf_edges: [plan_meta, plan_meta.settings.maf_edges]
        }

    PLINK_PREPARE_LDMS_PLAN_GCTA(
        ch_plan_requests.genotypes,
        ch_plan_requests.region_kb,
        ch_plan_requests.ld_bins,
        ch_plan_requests.maf_edges,
    )

    def ch_plans_by_key = PLINK_PREPARE_LDMS_PLAN_GCTA.out.plan.map { plan_meta, strata_manifest, snp_group_files -> [plan_meta.key, strata_manifest, snp_group_files] }

    // LDMS remains one base-family artifact. Its reusable product is the ordered GRM family; each native
    // consumer writes its own task-local MGRM control list from the explicit prefix order.
    def ch_ldms_builds = ch_base_by_type.gcta_ldms
        .map { _key, matrix_meta, base, _weights_file, parts -> [base.plan_key, base.view_compatibility_key, matrix_meta, parts] }
        .combine(ch_plans_by_key, by: 0)
        .map { _plan_key, view_key, matrix_meta, parts, strata_manifest, snp_group_files -> [view_key, matrix_meta, parts, strata_manifest, snp_group_files] }
        .combine(ch_plink1_by_view_key, by: 0)
        .map { _view_key, matrix_meta, parts, strata_manifest, snp_group_files, bed, bim, fam -> [matrix_meta, parts, strata_manifest, snp_group_files, bed, bim, fam] }

    def ch_ldms_manifests = ch_ldms_builds
        .map { matrix_meta, _parts, _strata_manifest, _snp_group_files, bed, _bim, _fam -> [matrix_meta, bed.baseName] }
        .collectFile { matrix_meta, stem -> ["${matrix_meta.id}.mbfile", "${stem}\n"] }
        .map { manifest -> [manifest.baseName, manifest] }

    def ch_ldms_inputs = ch_ldms_builds
        .map { matrix_meta, parts, strata_manifest, snp_group_files, bed, bim, fam -> [matrix_meta.id, matrix_meta, parts, strata_manifest, snp_group_files, bed, bim, fam] }
        .join(ch_ldms_manifests, failOnDuplicate: true, failOnMismatch: true)
        .multiMap { _matrix_id, matrix_meta, parts, strata_manifest, snp_group_files, bed, bim, fam, manifest ->
            plan: [matrix_meta, strata_manifest, snp_group_files]
            genotypes: [matrix_meta, manifest, bed, bim, fam]
            n_parts: [matrix_meta, parts]
        }

    PLINK_PREPARE_GRM_LDMS_GCTA(
        ch_ldms_inputs.plan,
        ch_ldms_inputs.genotypes,
        ch_ldms_inputs.n_parts,
    )

    // MPH matrices are built from the PLINK 1 bundle directly rather than from the prepared PGEN view, and
    // never from a GCTA matrix: the two layouts are mutually unreadable and neither tool rejects the other's.
    def ch_mph_dense_builds = ch_base_by_type.mph_dense
        .map { _key, matrix_meta, base, _weights_file, _parts -> [base.view_compatibility_key, matrix_meta] }
        .combine(ch_plink1_by_view_key, by: 0)
        .map { _view_key, matrix_meta, bed, bim, fam -> [matrix_meta, bed, bim, fam] }

    PLINK_PREPARE_GRM_MPH(ch_mph_dense_builds, autosome_count)

    def ch_mph_ldms_builds = ch_base_by_type.mph_ldms
        .map { _key, matrix_meta, base, _weights_file, _parts -> [base.plan_key, base.view_compatibility_key, matrix_meta] }
        .combine(ch_plans_by_key, by: 0)
        .map { _plan_key, view_key, matrix_meta, strata_manifest, snp_group_files -> [view_key, matrix_meta, strata_manifest, snp_group_files] }
        .combine(ch_plink1_by_view_key, by: 0)
        .multiMap { _view_key, matrix_meta, strata_manifest, snp_group_files, bed, bim, fam ->
            plan: [matrix_meta, strata_manifest, snp_group_files]
            genotypes: [matrix_meta, bed, bim, fam]
        }

    PLINK_PREPARE_GRM_LDMS_MPH(ch_mph_ldms_builds.plan, ch_mph_ldms_builds.genotypes, autosome_count)

    // CALCKINS sees only LDAK base requests. Optional filtering/subsetting is supplied as a child request whose
    // parent key selects the completed base, so unrestricted and filtered consumers share construction.
    def ch_ldak_inputs = ch_base_by_type.ldak_kinship
        .map { key, matrix_meta, base, weights_file, _parts -> [base.view_compatibility_key, key, matrix_meta, base, weights_file] }
        .combine(ch_plink1_by_view_key, by: 0)
        .map { _view_key, key, matrix_meta, base, weights_file, bed, bim, fam -> [key, matrix_meta, base, weights_file, bed, bim, fam] }
        .multiMap { _key, matrix_meta, base, weights_file, bed, bim, fam ->
            genotypes: [matrix_meta, bed, bim, fam, base.settings.power]
            weights: [matrix_meta, weights_file]
        }

    def ch_ldak_subset_requests = ch_requests
        .filter { _selected_key, _meta, request, _weights_file -> request.derived && request.derived.type == 'ldak_subset' }
        .map { _selected_key, _meta, request, _weights_file -> [request.derived.parent_key, buildRelatednessArtifactMeta(request.derived)] }
        .unique { _parent_key, derived_meta -> derived_meta.key }

    PLINK_PREPARE_GRM_LDAK(
        ch_ldak_inputs.genotypes,
        ch_ldak_inputs.weights,
        ch_ldak_subset_requests,
    )

    // Fan every selected final artifact back to its unchanged focal analysis/request metadata.
    def ch_gcta_dense = ch_requests
        .filter { _selected_key, _meta, request, _weights_file -> request.kind == 'gcta_dense' }
        .map { selected_key, meta, _request, _weights_file -> [selected_key, meta.relationship_id ? meta + [matrix_key: selected_key] : meta] }
        .combine(
            PLINK_PREPARE_GRM_GCTA.out.grm_files.map { matrix_meta, grm_files -> [matrix_meta.key, grm_files] },
            by: 0
        )
        .map { _selected_key, meta, grm_files -> [meta, grm_files] }

    def ch_gcta_sparse = ch_requests
        .filter { _selected_key, _meta, request, _weights_file -> request.kind == 'gcta_sparse' }
        .map { selected_key, meta, _request, _weights_file -> [selected_key, meta] }
        .combine(
            GCTA_MAKEBKSPARSE.out.sparse_grm_files.map { matrix_meta, sparse_grm_files -> [matrix_meta.key, sparse_grm_files] },
            by: 0
        )
        .map { _selected_key, meta, sparse_grm_files -> [meta, sparse_grm_files] }

    def ch_gcta_ldms = ch_requests
        .filter { _selected_key, _meta, request, _weights_file -> request.kind == 'gcta_ldms' }
        .map { selected_key, meta, _request, _weights_file -> [selected_key, meta.relationship_id ? meta + [matrix_key: selected_key] : meta] }
        .combine(
            PLINK_PREPARE_GRM_LDMS_GCTA.out.grm_family.map { matrix_meta, grm_files, grm_prefixes -> [matrix_meta.key, grm_files, grm_prefixes] },
            by: 0
        )
        .map { _selected_key, meta, grm_files, grm_prefixes -> [meta, grm_files, grm_prefixes] }

    // The MPH families carry their identity in a tuple position rather than in analysis metadata. The matrix
    // identity includes the plan and declared predictor counts read from the completed construction outputs.
    def ch_mph_dense_artifacts = PLINK_PREPARE_GRM_MPH.out.grm_files
        .map { matrix_meta, grm_files -> [matrix_meta.key, matrix_meta, grm_files] }
        .join(
            PLINK_PREPARE_GRM_MPH.out.counts.map { matrix_meta, counts -> [matrix_meta.key, counts] },
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { key, matrix_meta, grm_files, counts -> [key, buildMphMatrixIdentity(matrix_meta, counts), grm_files] }

    def ch_mph_dense = ch_requests
        .filter { _selected_key, _meta, request, _weights_file -> request.kind == 'mph_dense' }
        .map { selected_key, meta, _request, _weights_file -> [selected_key, meta] }
        .combine(ch_mph_dense_artifacts, by: 0)
        .map { _selected_key, meta, matrix_identity, grm_files -> [meta, matrix_identity, grm_files] }

    def ch_mph_ldms_artifacts = PLINK_PREPARE_GRM_LDMS_MPH.out.grm_family
        .map { matrix_meta, grm_files, grm_prefixes -> [matrix_meta.key, matrix_meta, grm_files, grm_prefixes] }
        .join(
            PLINK_PREPARE_GRM_LDMS_MPH.out.counts.map { matrix_meta, counts -> [matrix_meta.key, counts] },
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { key, matrix_meta, grm_files, grm_prefixes, counts -> [key, buildMphMatrixIdentity(matrix_meta, counts), grm_files, grm_prefixes] }

    def ch_mph_ldms = ch_requests
        .filter { _selected_key, _meta, request, _weights_file -> request.kind == 'mph_ldms' }
        .map { selected_key, meta, _request, _weights_file -> [selected_key, meta] }
        .combine(ch_mph_ldms_artifacts, by: 0)
        .map { _selected_key, meta, matrix_identity, grm_files, grm_prefixes -> [meta, matrix_identity, grm_files, grm_prefixes] }

    def ch_ldak_unfiltered = ch_requests
        .filter { _selected_key, _meta, request, _weights_file -> request.kind == 'ldak_kinship' && !request.derived }
        .map { selected_key, meta, _request, _weights_file -> [selected_key, meta] }
        .combine(
            PLINK_PREPARE_GRM_LDAK.out.base_grm.map { matrix_meta, grm_files -> [matrix_meta.key, grm_files] },
            by: 0
        )
        .map { selected_key, meta, grm_files -> [meta, selected_key, grm_files, []] }

    def ch_ldak_filtered = ch_requests
        .filter { _selected_key, _meta, request, _weights_file -> request.kind == 'ldak_kinship' && request.derived }
        .map { selected_key, meta, _request, _weights_file -> [selected_key, meta] }
        .combine(
            PLINK_PREPARE_GRM_LDAK.out.subset_grm.map { matrix_meta, grm_files -> [matrix_meta.key, grm_files] },
            by: 0
        )
        .combine(
            PLINK_PREPARE_GRM_LDAK.out.filtered_list.map { matrix_meta, keep, _lose -> [matrix_meta.key, keep] },
            by: 0
        )
        .map { selected_key, meta, grm_files, keep -> [meta, selected_key, grm_files, keep] }

    def ch_ldak_kinship = ch_ldak_unfiltered.mix(ch_ldak_filtered)

    // The component plan is published as its own artifact family, addressed by its own key, because it is now
    // shared by two matrix families and belongs to neither.
    def ch_ldms_plan_artifacts = PLINK_PREPARE_LDMS_PLAN_GCTA.out.plan
        .map { plan_meta, strata_manifest, snp_group_files -> [plan_meta.key, plan_meta, strata_manifest, snp_group_files] }
        .join(
            PLINK_PREPARE_LDMS_PLAN_GCTA.out.ld_scores.map { plan_meta, ld_scores -> [plan_meta.key, ld_scores] },
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { _key, plan_meta, strata_manifest, snp_group_files, ld_scores -> [plan_meta, ld_scores, strata_manifest, snp_group_files] }

    emit:
    gcta_dense          = ch_gcta_dense // channel: [ val(meta), path(grm_files) ]
    gcta_sparse         = ch_gcta_sparse // channel: [ val(meta), path(sparse_grm_files) ]
    gcta_ldms           = ch_gcta_ldms // channel: [ val(meta), path(grm_files), val(grm_prefixes) ]
    gcta_ldms_artifacts = PLINK_PREPARE_GRM_LDMS_GCTA.out.grm_family // channel: [ val(matrix_meta), path(grm_files), val(grm_prefixes) ], one per base key
    ldak_kinship        = ch_ldak_kinship // channel: [ val(meta), val(artifact_key), path(grm_files), path(keep) ]
    mph_dense           = ch_mph_dense // channel: [ val(meta), val(matrix_identity), path(grm_files) ]
    mph_ldms            = ch_mph_ldms // channel: [ val(meta), val(matrix_identity), path(grm_files), val(grm_prefixes) ]
    mph_ldms_artifacts  = PLINK_PREPARE_GRM_LDMS_MPH.out.grm_family // channel: [ val(matrix_meta), path(grm_files), val(grm_prefixes) ], one per base key
    ldms_plan_artifacts = ch_ldms_plan_artifacts // channel: [ val(plan_meta), path(ld_scores), path(strata_manifest), [ path(snp_group_file), ... ] ], one per plan key
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Choose one cohort's bundle members when several cohorts share one genotype view. The cohort id is the only
// stable ordering available — the members are byte-identical and their paths are what must not decide — so
// the lowest one wins, which makes the choice a function of the manifest rather than of task completion
// order. `members` is the per-member list of lists in tuple order; the answer is one member per list.
def selectLowestCohort(cohort_ids, members) {
    def chosen = cohort_ids.indexOf(cohort_ids.min())
    return members.collect { member -> member[chosen] }
}

def canonicaliseDeclaredValue(value) {
    return canonicaliseScientificValue(value)
}

def canonicaliseIdentifier(value) {
    return canonicaliseScientificIdentifier(value)
}

def buildRelatednessMatrixKey(identity, settings) {
    return buildScientificArtifactKey(identity, settings)
}

def buildRelatednessArtifactMeta(artifact) {
    def meta = [
        id: "${artifact.type}.${artifact.key}",
        key: artifact.key,
        layer: artifact.layer,
        artifact_type: artifact.type,
        settings: artifact.settings,
    ]
    if (artifact.plan) {
        meta = meta + [plan_key: artifact.plan.key, plan_settings: artifact.plan.settings]
    }
    return artifact.layer == 'base'
        ? meta + [view_compatibility_key: artifact.view_compatibility_key]
        : meta + [parent_key: artifact.parent_key]
}

def buildLdmsPlanMeta(plan) {
    return [
        id: "${plan.type}.${plan.key}",
        key: plan.key,
        layer: plan.layer,
        artifact_type: plan.type,
        settings: plan.settings,
        view_compatibility_key: plan.view_compatibility_key,
    ]
}

// Read the declared predictor counts an MPH builder emitted beside its matrices from the SNP-information
// file at build time. These are distinct from MPH's post-quality-control counts in the fit's result. The `.grm.bin` header's
// `sum2pq` is deliberately not used for either: it is a `float32` sum of 2pq weights that coincides with a
// variant count only for the 0/1 columns this pipeline happens to write, and it is not exact above 2^24.
def readMphPredictorCounts(counts_file) {
    def reserved = ['__unassigned__', '__bim_rows__']
    def records = counts_file
        .readLines()
        .findAll { row -> row.trim() }
        .drop(1)
        .collect { row -> row.split('\t', -1).toList() }
    def component_records = records.findAll { record -> !(record[0] in reserved) }
    def totals = records.findAll { record -> record[0] in reserved }.collectEntries { record -> [(record[0]): record[1] as int] }
    def components = component_records
        .withIndex()
        .collect { record, index -> [ordinal: index + 1, name: record[0], plan_predictor_count: record[1] as int] }
    return [
        components: components,
        variant_counts: [
            bim_rows: totals['__bim_rows__'],
            assigned: components.sum { component -> component.plan_predictor_count } ?: 0,
            unassigned: totals['__unassigned__'],
        ],
    ]
}

def buildMphMatrixIdentity(matrix_meta, counts_file) {
    def counts = readMphPredictorCounts(counts_file)
    return [
        key: matrix_meta.key,
        layer: matrix_meta.layer,
        artifact_type: matrix_meta.artifact_type,
        settings: matrix_meta.settings,
        plan_key: matrix_meta.plan_key,
        plan_settings: matrix_meta.plan_settings,
        components: counts.components,
        variant_counts: counts.variant_counts,
    ]
}

def getRelatednessMatrixKinds(meta) {
    def capabilities = getMethodCapabilities()
    def selected = meta.relationship_id
        ? [meta.method] as Set
        : ((meta.association_methods ?: []) + (meta.heritability_methods ?: [])) as Set
    return getMatrixKindContract()
        .keySet()
        .toList()
        .findAll { kind ->
            selected.any { method -> capabilities[method] && capabilities[method].matrix_kind == kind }
        }
}

def getLdakWeightsIdentity(weights_file, weights_policy = 'equal') {
    if (!weights_file) {
        return [mode: weights_policy]
    }
    return [mode: 'provided', sha256: digestFileBytes(weights_file)]
}

def getMethodResourceIdentity(resource) {
    if (!resource) {
        return [mode: 'all']
    }
    return [mode: 'file', sha256: digestFileBytes(resource)]
}

// The three settings that define one LD-by-MAF component plan. They stay in the `gcta` option family because
// GCTA's own LD-score pass builds the plan, and both stratified matrix families read them from here, so a GCTA
// and an MPH request that declare the same three settings share one plan rather than computing it twice.
def getLdmsPlanSettings(meta) {
    return meta.relationship_id
        ? meta.matrix_settings
        : [
            ld_score_region_kb: meta.method_options.gcta.ld_score_region_kb,
            ld_bins: meta.method_options.gcta.ld_bins,
            maf_edges: meta.method_options.gcta.ldms_maf_edges,
        ]
}

def getRelatednessBaseSettings(meta, requested_kind, weights_identity = [mode: 'equal'], gcta_extract_identity = [mode: 'all']) {
    if (requested_kind in ['gcta_dense', 'gcta_sparse']) {
        if (meta.relationship_id) {
            return meta.matrix_settings
        }
        def settings = [:]
        if (meta.method_options.gcta.grm_maf != null) {
            settings.maf = meta.method_options.gcta.grm_maf
        }
        if (gcta_extract_identity.mode == 'file') {
            settings.extract = gcta_extract_identity
        }
        return settings
    }
    if (requested_kind == 'gcta_ldms') {
        return getLdmsPlanSettings(meta)
    }
    if (requested_kind == 'ldak_kinship') {
        return [
            model: meta.method_options.ldak.model,
            power: meta.method_options.ldak.power,
            weights: weights_identity,
        ]
    }
    // The MPH families name their on-disk format rather than a construction option, because MPH's matrix
    // takes none: no MAF filter, no extract list and no scaling exponent enter its construction. The dense
    // family additionally names its variant universe symbolically, so a later change of universe cannot reuse
    // today's bytes under today's key. The stratified family says nothing about the universe here because its
    // plan is autosome-only by construction and the plan key already carries the settings that fixed it.
    if (requested_kind == 'mph_dense') {
        return [format: 'mph_grm', variant_universe: 'autosomes']
    }
    if (requested_kind == 'mph_ldms') {
        return [format: 'mph_grm']
    }
}

def getRelatednessDerivedSettings(meta, requested_kind) {
    if (requested_kind == 'gcta_sparse') {
        return [cutoff: meta.method_options.gcta.sparse_cutoff]
    }
    if (requested_kind == 'ldak_kinship' && meta.method_options.ldak.relatedness_filter) {
        return [relatedness_filter: true]
    }
    return [:]
}

// The view-compatibility key arrives as an argument rather than being derived here: which of a cohort's two
// views a request is compatible with is decided by the kind's declared genotype bundle, and only the caller
// holds both key channels. Nothing in this function reads a staged genotype name.
def buildRelatednessArtifactRequest(meta, view_compatibility_key, kind, weights_identity = [mode: 'equal'], gcta_extract_identity = [mode: 'all'], gcta_extract = []) {
    def base_type = kind == 'gcta_sparse' ? 'gcta_dense' : kind
    def base_settings = getRelatednessBaseSettings(meta, kind, weights_identity, gcta_extract_identity)

    // A stratified base depends on the plan, so it names the plan key in its identity rather than repeating
    // the plan settings. `gcta_ldms` is the exception on purpose: its key was published before the plan became
    // its own artifact, and it keeps the original formula so no published matrix directory moves.
    def plan = null
    if (base_type in ['gcta_ldms', 'mph_ldms']) {
        def plan_settings = getLdmsPlanSettings(meta)
        plan = [
            layer: 'plan',
            type: 'ldms_component_plan',
            key: buildScientificArtifactKey(
                [layer: 'plan', type: 'ldms_component_plan', view_compatibility_key: view_compatibility_key],
                plan_settings,
            ),
            settings: plan_settings,
            view_compatibility_key: view_compatibility_key,
        ]
    }

    def base_key = base_type == 'mph_ldms'
        ? buildScientificArtifactKey([layer: 'base', type: base_type, plan_key: plan.key], base_settings)
        : buildScientificArtifactKey([layer: 'base', type: base_type, view_compatibility_key: view_compatibility_key], base_settings)
    def base = [
        layer: 'base',
        type: base_type,
        key: base_key,
        view_compatibility_key: view_compatibility_key,
        settings: base_settings,
    ]
    if (plan) {
        base.plan = plan
        base.plan_key = plan.key
    }
    if (base_type == 'gcta_dense') {
        base.gcta_extract = gcta_extract
    }

    def derived = [:]
    if (kind == 'gcta_sparse') {
        def settings = getRelatednessDerivedSettings(meta, kind)
        derived = [
            layer: 'derived',
            type: 'gcta_sparse',
            parent_key: base.key,
            settings: settings,
            key: buildScientificArtifactKey([layer: 'derived', type: 'gcta_sparse', parent_key: base.key], settings),
        ]
    }
    if (kind == 'ldak_kinship' && meta.method_options.ldak.relatedness_filter) {
        def settings = getRelatednessDerivedSettings(meta, kind)
        derived = [
            layer: 'derived',
            type: 'ldak_subset',
            parent_key: base.key,
            settings: settings,
            key: buildScientificArtifactKey([layer: 'derived', type: 'ldak_subset', parent_key: base.key], settings),
        ]
    }

    return [
        kind: kind,
        base: base,
        derived: derived,
        selected_key: derived ? derived.key : base.key,
    ]
}
