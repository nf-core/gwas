// Build each distinct relatedness matrix once and fan it out to every analysis unit that needs it.
// Every component reports on the run-wide versions topic, so this subworkflow emits no versions.

// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
include { PLINK_PREPARE_GRM_GCTA      } from '../plink_prepare_grm_gcta/main'
include { PLINK_PREPARE_GRM_LDMS_GCTA } from '../plink_prepare_grm_ldms_gcta/main'
include { PLINK_PREPARE_GRM_LDAK      } from '../plink_prepare_grm_ldak/main'

// MODULE: Local to the pipeline
include { GCTA_MAKEBKSPARSE           } from '../../../modules/local/gcta/makebksparse/main'

// FUNCTION: Local to the pipeline
include { digestFileBytes             } from '../utils_nfcore_gwas_pipeline'
include { digestIdentityText          } from '../utils_nfcore_gwas_pipeline'
include { getMethodCapabilities       } from '../validate_gwas_input'

workflow PREPARE_RELATEDNESS_MATRICES {
    take:
    ch_analyses // channel: [ val(meta), [ path(genotype_file), ... ], path(ldak_weights) ], [] when absent
    ch_cohort_genotypes // channel: [ val(cohort_meta), path(pgen), path(psam), path(pvar) ]
    ch_plink1_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam) ], analyses needing PLINK 1
    gcta_grm_parts // channel: val(gcta_grm_parts), run/profile GCTA GRM partition count

    main:

    // Relatedness matrices are the most expensive artefacts built here. Hoisting construction above the
    // heritability routes lets analyses with the same cohort and kinship model reuse one matrix.
    //
    // One request per analysis unit per matrix it needs, keyed by the reuse digest. An analysis selecting
    // three estimators that share a matrix contributes three requests carrying one key; an analysis
    // selecting none contributes nothing.
    //
    def ch_requests = ch_analyses.flatMap { meta, genotype_files, ldak_weights ->
        getRelatednessMatrixKinds(meta).collect { kind ->
            def weights_policy = meta.method_options.ldak.weights_policy
            def weights_identity = kind == 'ldak_kinship' ? getLdakWeightsIdentity(ldak_weights, weights_policy) : [mode: 'equal']
            def gcta_extract = kind == 'gcta_dense' ? meta.method_options.gcta.grm_extract : []
            def extract_identity = kind == 'gcta_dense' ? getMethodResourceIdentity(gcta_extract) : [mode: 'all']
            def request = buildRelatednessMatrixRequest(meta, genotype_files, kind, weights_identity, extract_identity)
            def weights_file = kind == 'ldak_kinship' ? ldak_weights ?: [] : []
            [request.key, meta, request, weights_file]
        }
    }

    //
    // GCTA partitioning is one run/profile value. It is deliberately not reconciled per matrix because no
    // analysis can declare it, and it remains outside every content identity because partitioning changes
    // task shape rather than merged matrix bytes.
    //

    //
    // One build per distinct key. Both sides of this join are reductions of `ch_requests`, so the key sets
    // are equal and unique by construction and the strict form is what says so if that ever stops being
    // true. The matrix identity is built only from key components: nothing analysis-derived reaches the
    // processes below, which is what makes per-trait missingness — and any other per-trait value —
    // structurally incapable of fragmenting the matrix.
    //
    def ch_matrices = ch_requests
        .map { key, _meta, request, weights_file -> [key, request, weights_file] }
        .unique { key, _request, _weights_file -> key }
        .map { key, request, weights_file ->
            def matrix_meta = [
                id: "${request.cohort}.${request.kind}.${key}",
                cohort: request.cohort,
                kind: request.kind,
                key: key,
                settings: request.settings,
            ]
            [key, matrix_meta, request, weights_file, gcta_grm_parts]
        }

    //
    // The cohort seam, again a combining operator and for the same reason it is one in
    // PREPARE_COHORT_GENOTYPES: one cohort carries one prepared bundle and may carry several matrices, so
    // `join` would emit only the first matrix of each cohort and its strict form would raise on the second.
    // Exactly one bundle exists per cohort key, so the product is exactly one element per matrix.
    //
    def ch_matrix_genotypes = ch_matrices
        .filter { _key, matrix_meta, _request, _weights_file, _parts ->
            matrix_meta.kind in ['gcta_dense', 'gcta_sparse', 'gcta_ldms']
        }
        .map { _key, matrix_meta, request, _weights_file, parts ->
            [matrix_meta.cohort, matrix_meta, parts, request.gcta_extract ?: []]
        }
        .combine(
            ch_cohort_genotypes.map { cohort_meta, pgen, psam, pvar -> [cohort_meta.id, pgen, psam, pvar] },
            by: 0
        )
        .map { _cohort_id, matrix_meta, parts, gcta_extract, pgen, psam, pvar ->
            [matrix_meta, parts, gcta_extract, pgen, psam, pvar]
        }

    def ch_by_kind = ch_matrix_genotypes.branch { matrix_meta, _parts, _gcta_extract, _pgen, _psam, _pvar ->
        gcta_dense: matrix_meta.kind == 'gcta_dense'
        gcta_sparse: matrix_meta.kind == 'gcta_sparse'
        gcta_ldms: matrix_meta.kind == 'gcta_ldms'
    }

    //
    // GCTA resolves every manifest entry from one prefix in the task working directory. Derive that
    // prefix from the staged primary PGEN file; matching PSAM/PVAR basenames are the native caller
    // contract, so invalid bundles fail in GCTA. `collectFile` keeps the one-line manifest cacheable.
    //
    def ch_dense_builds = ch_by_kind.gcta_dense.mix(ch_by_kind.gcta_sparse)

    def ch_dense_manifests = ch_dense_builds
        .map { matrix_meta, _parts, _gcta_extract, pgen, _psam, _pvar -> [matrix_meta, pgen.baseName] }
        .collectFile { matrix_meta, stem -> ["${matrix_meta.id}.mpfile", "${stem}\n"] }
        .map { manifest -> [manifest.baseName, manifest] }

    def ch_dense_inputs = ch_dense_builds
        .map { matrix_meta, parts, gcta_extract, pgen, psam, pvar -> [matrix_meta.id, matrix_meta, parts, gcta_extract, pgen, psam, pvar] }
        .join(ch_dense_manifests, failOnDuplicate: true, failOnMismatch: true)
        .multiMap { _matrix_id, matrix_meta, parts, gcta_extract, pgen, psam, pvar, manifest ->
            genotypes: [matrix_meta, manifest, pgen, pvar, psam]
            snp_group: [matrix_meta, gcta_extract]
            n_parts: [matrix_meta, parts]
        }

    //
    // SUBWORKFLOW: Build the dense GCTA relatedness matrix in partitions and merge them
    //
    PLINK_PREPARE_GRM_GCTA(
        ch_dense_inputs.genotypes,
        ch_dense_inputs.snp_group,
        ch_dense_inputs.n_parts,
    )

    //
    // MODULE: Convert the dedicated dense base of each fastGWA request to GCTA's sparse format.
    //
    // The sparse cutoff lives inside the matrix identity and reaches the component as an explicit input.
    // It is not hidden in ext.args, so changing it necessarily changes both the reuse key and the command.
    //
    def ch_sparse_inputs = PLINK_PREPARE_GRM_GCTA.out.grm_files
        .filter { matrix_meta, _grm_files -> matrix_meta.kind == 'gcta_sparse' }
        .multiMap { matrix_meta, grm_files ->
            grm: [matrix_meta, grm_files]
            cutoff: matrix_meta.settings.cutoff
        }

    GCTA_MAKEBKSPARSE(
        ch_sparse_inputs.grm,
        ch_sparse_inputs.cutoff,
    )

    //
    // SUBWORKFLOW: Build one LD-by-MAF-stratified matrix family per distinct LDMS reuse key.
    //
    // LD-score calculation consumes PLINK 1, so the matrix request is paired with the lazy cohort
    // derivative prepared upstream. The derivative is emitted per requesting analysis; reducing it to
    // one record per cohort here preserves one matrix build per key without a second conversion.
    //
    def ch_plink1_cohorts = ch_plink1_genotypes
        .map { meta, bed, bim, fam -> [meta.cohort, bed, bim, fam] }
        .unique { cohort, _bed, _bim, _fam -> cohort }

    def ch_ldms_builds = ch_by_kind.gcta_ldms
        .map { matrix_meta, parts, _gcta_extract, _pgen, _psam, _pvar -> [matrix_meta.cohort, matrix_meta, parts] }
        .combine(ch_plink1_cohorts, by: 0)
        .map { _cohort, matrix_meta, parts, bed, bim, fam -> [matrix_meta, parts, bed, bim, fam] }

    // Derive the LDMS PLINK 1 prefix from the staged primary BED file under the same native bundle contract.
    def ch_ldms_manifests = ch_ldms_builds
        .map { matrix_meta, _parts, bed, _bim, _fam -> [matrix_meta, bed.baseName] }
        .collectFile { matrix_meta, stem -> ["${matrix_meta.id}.mbfile", "${stem}\n"] }
        .map { manifest -> [manifest.baseName, manifest] }

    def ch_ldms_inputs = ch_ldms_builds
        .map { matrix_meta, parts, bed, bim, fam -> [matrix_meta.id, matrix_meta, parts, bed, bim, fam] }
        .join(ch_ldms_manifests, failOnDuplicate: true, failOnMismatch: true)
        .multiMap { _matrix_id, matrix_meta, parts, bed, bim, fam, manifest ->
            genotypes: [matrix_meta, manifest, bed, bim, fam]
            ld_score_region_kb: [matrix_meta, matrix_meta.settings.ld_score_region_kb]
            ld_bins: [matrix_meta, matrix_meta.settings.ld_bins]
            maf_edges: [matrix_meta, matrix_meta.settings.maf_edges]
            n_parts: [matrix_meta, parts]
        }

    PLINK_PREPARE_GRM_LDMS_GCTA(
        ch_ldms_inputs.genotypes,
        ch_ldms_inputs.ld_score_region_kb,
        ch_ldms_inputs.ld_bins,
        ch_ldms_inputs.maf_edges,
        ch_ldms_inputs.n_parts,
    )

    //
    // SUBWORKFLOW: Build one LDAK kinship matrix per distinct scientific identity. Relatedness filtering is
    // part of that identity, so an unrestricted request and its unrelated-subset counterpart build and
    // publish independently while analyses agreeing on the filter still share construction.

    // PREPARE_COHORT_GENOTYPES already built one PLINK 1 derivative per requesting cohort and fanned it out
    // with focal analysis metadata. Collapse those identical fan-out records back to one cohort bundle before
    // combining them with the possibly-many matrix keys of that cohort.
    def ch_ldak_inputs = ch_matrices
        .filter { _key, matrix_meta, _request, _weights_file, _parts -> matrix_meta.kind == 'ldak_kinship' }
        .map { key, matrix_meta, request, weights_file, _parts -> [matrix_meta.cohort, key, matrix_meta, request, weights_file] }
        .combine(ch_plink1_cohorts, by: 0)
        .map { _cohort, key, matrix_meta, request, weights_file, bed, bim, fam -> [key, matrix_meta, request, weights_file, bed, bim, fam] }
        .multiMap { _key, matrix_meta, request, weights_file, bed, bim, fam ->
            genotypes: [matrix_meta, bed, bim, fam, request.settings.power]
            weights: [matrix_meta, weights_file]
            filter_relatedness: [matrix_meta, request.filter_relatedness]
        }

    PLINK_PREPARE_GRM_LDAK(
        ch_ldak_inputs.genotypes,
        ch_ldak_inputs.weights,
        ch_ldak_inputs.filter_relatedness,
    )

    //
    // The key-to-analysis seam. One key to many analysis units, so a combining operator again, and because
    // exactly one matrix exists per key the product is exactly one element per requesting analysis unit,
    // carrying the focal analysis meta rather than the matrix meta.
    //
    def ch_gcta_dense = ch_requests
        .filter { _key, _meta, request, _weights_file -> request.kind == 'gcta_dense' }
        .map { key, meta, _request, _weights_file -> [key, meta] }
        .combine(
            PLINK_PREPARE_GRM_GCTA.out.grm_files.map { matrix_meta, grm_files -> [matrix_meta.key, grm_files] },
            by: 0
        )
        .map { _key, meta, grm_files -> [meta, grm_files] }

    def ch_gcta_sparse = ch_requests
        .filter { _key, _meta, request, _weights_file -> request.kind == 'gcta_sparse' }
        .map { key, meta, _request, _weights_file -> [key, meta] }
        .combine(
            GCTA_MAKEBKSPARSE.out.sparse_grm_files.map { matrix_meta, sparse_grm_files -> [matrix_meta.key, sparse_grm_files] },
            by: 0
        )
        .map { _key, meta, sparse_grm_files -> [meta, sparse_grm_files] }

    def ch_gcta_ldms = ch_requests
        .filter { _key, _meta, request, _weights_file -> request.kind == 'gcta_ldms' }
        .map { key, meta, _request, _weights_file -> [key, meta] }
        .combine(
            PLINK_PREPARE_GRM_LDMS_GCTA.out.mgrm_bundle.map { matrix_meta, mgrm, grm_files -> [matrix_meta.key, mgrm, grm_files] },
            by: 0
        )
        .map { _key, meta, mgrm, grm_files -> [meta, mgrm, grm_files] }

    // Unrestricted analyses take the all-sample emission. Restricted analyses take the derived unrelated
    // matrix and carry its keep list onward to the estimator; their matrix keys remain distinct because the
    // requested filtering policy changes the scientific matrix identity.
    def ch_ldak_unfiltered = ch_requests
        .filter { _key, _meta, request, _weights_file -> request.kind == 'ldak_kinship' && !request.filter_relatedness }
        .map { key, meta, _request, _weights_file -> [key, meta] }
        .combine(
            PLINK_PREPARE_GRM_LDAK.out.unfiltered_grm.map { matrix_meta, grm_files -> [matrix_meta.key, grm_files] },
            by: 0
        )
        .map { _key, meta, grm_files -> [meta, grm_files, []] }

    def ch_ldak_filtered = ch_requests
        .filter { _key, _meta, request, _weights_file -> request.kind == 'ldak_kinship' && request.filter_relatedness }
        .map { key, meta, _request, _weights_file -> [key, meta] }
        .combine(
            PLINK_PREPARE_GRM_LDAK.out.analysis_grm.map { matrix_meta, grm_files -> [matrix_meta.key, grm_files] },
            by: 0
        )
        .combine(
            PLINK_PREPARE_GRM_LDAK.out.filtered_list.map { matrix_meta, keep, _lose -> [matrix_meta.key, keep] },
            by: 0
        )
        .map { _key, meta, grm_files, keep -> [meta, grm_files, keep] }

    def ch_ldak_kinship = ch_ldak_unfiltered.mix(ch_ldak_filtered)

    emit:
    gcta_dense   = ch_gcta_dense // channel: [ val(meta), path(grm_files) ]
    gcta_sparse  = ch_gcta_sparse // channel: [ val(meta), path(sparse_grm_files) ]
    gcta_ldms    = ch_gcta_ldms // channel: [ val(meta), path(mgrm), path(grm_files) ]
    ldak_kinship = ch_ldak_kinship // channel: [ val(meta), path(grm_files), path(keep) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Declared scientific values use numeric canonicalisation so equal quantities share one matrix.
def canonicaliseDeclaredValue(value) {
    if (value == null) {
        return 'null'
    }
    if (value instanceof Boolean) {
        return value ? 'true' : 'false'
    }
    if (value instanceof Map) {
        return '{' + value.sort { entry -> entry.key }.collect { name, entry -> "${name}=${canonicaliseDeclaredValue(entry)}" }.join(',') + '}'
    }
    if (value instanceof Collection) {
        return '[' + value.collect { entry -> canonicaliseDeclaredValue(entry) }.join(',') + ']'
    }
    if (value instanceof Number) {
        return new BigDecimal(value.toString()).stripTrailingZeros().toPlainString()
    }
    def text = value.toString().trim()
    return text ==~ /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/
        ? new BigDecimal(text).stripTrailingZeros().toPlainString()
        : text
}

// Identifiers are names, not quantities: cohorts `01` and `1` must never collapse onto one matrix.
def canonicaliseIdentifier(value) {
    if (value == null) {
        return 'null'
    }
    if (value instanceof Map) {
        return '{' + value.sort { entry -> entry.key }.collect { name, entry -> "${name}=${canonicaliseIdentifier(entry)}" }.join(',') + '}'
    }
    if (value instanceof Collection) {
        return '[' + value.collect { entry -> canonicaliseIdentifier(entry) }.join(',') + ']'
    }
    return value.toString().trim()
}

// Relatedness keys retain their exact identity/settings serialisation and delegate only the SHA primitive.
def buildRelatednessMatrixKey(identity, settings) {
    def rendered = identity.collectEntries { name, value -> [(name): canonicaliseIdentifier(value)] } + [settings: canonicaliseDeclaredValue(settings)]
    def canonical = rendered
        .sort { entry -> entry.key }
        .collect { name, text -> "${name}=${text}" }
        .join('\n')
    return digestIdentityText(canonical)
}

def getRelatednessMatrixKinds(meta) {
    def capabilities = getMethodCapabilities()
    def selected = ((meta.association_methods ?: []) + (meta.heritability_methods ?: [])) as Set
    return ['gcta_dense', 'gcta_ldms', 'gcta_sparse', 'ldak_kinship'].findAll { kind ->
        selected.any { method -> capabilities[method] && capabilities[method].matrix_kind == kind }
    }
}

// Absence policy remains relatedness-specific; present resources share the common byte digest.
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

// Matrix settings contain every scientific construction input and exclude estimator/execution controls.
def getRelatednessMatrixSettings(meta, kind, weights_identity = [mode: 'equal'], gcta_extract_identity = [mode: 'all']) {
    def method_options = meta.method_options
    def ldak_options = method_options.ldak
    if (kind == 'gcta_dense') {
        def settings = [:]
        if (method_options.gcta.grm_maf != null) {
            settings.maf = method_options.gcta.grm_maf
        }
        if (gcta_extract_identity.mode == 'file') {
            settings.extract = gcta_extract_identity
        }
        return settings
    }
    if (kind == 'gcta_ldms') {
        return [
            ld_score_region_kb: method_options.gcta.ld_score_region_kb,
            ld_bins: method_options.gcta.ld_bins,
            maf_edges: method_options.gcta.ldms_maf_edges,
        ]
    }
    if (kind == 'gcta_sparse') {
        return [cutoff: method_options.gcta.sparse_cutoff]
    }
    if (kind == 'ldak_kinship') {
        def settings = [
            model: ldak_options.model,
            power: ldak_options.power,
            weights: weights_identity,
        ]
        if (ldak_options.relatedness_filter) {
            settings.relatedness_filter = true
        }
        return settings
    }
    error("[nf-core/gwas] ERROR: no relatedness matrix settings are registered for kind '${kind}' requested by analysis unit '${meta.id}'")
}

def buildRelatednessMatrixRequest(meta, genotype_files, kind, weights_identity = [mode: 'equal'], gcta_extract_identity = [mode: 'all']) {
    def method_options = meta.method_options
    def settings = getRelatednessMatrixSettings(meta, kind, weights_identity, gcta_extract_identity)
    def identity = [
        cohort: meta.cohort,
        genotype_format: meta.genotype_format,
        genotypes: genotype_files.collect { genotype_file -> genotype_file.name }.sort(),
        kind: kind,
    ]
    def request = [
        kind: kind,
        cohort: meta.cohort,
        settings: settings,
        key: buildRelatednessMatrixKey(identity, settings),
    ]
    if (kind == 'gcta_dense') {
        request.gcta_extract = method_options.gcta.grm_extract
    }
    if (kind == 'ldak_kinship') {
        request.filter_relatedness = method_options.ldak.relatedness_filter
    }
    return request
}
