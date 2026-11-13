// The closed vocabularies every registry entry is written against. Routing and validation ask the registry
// which methods hold a capability instead of repeating a method-name list that silently goes stale when a
// method is added, so a new entry that omits or misspells a capability must fail the registry contract test
// rather than quietly fall out of a hardcoded list.
//
// `input_backend` names the representation the estimator genuinely consumes. A direct-genotype estimator is a
// distinct backend and must never be given a fabricated matrix kind to make it look like a GRM route.
// `sparse_grm` is a real backend rather than a variant of `dense_grm`: fastGWA streams genotypes against a
// sparse relatedness matrix. The direct backends are split by genotype representation, not by tool:
// `direct_plink1_genotypes` names an executable that reads BED/BIM/FAM only, `direct_plink_genotypes` one that
// selects the native flag from the staged primary extension.
//
// `mph_grm` and `mph_grm_family` are separate backends from `dense_grm` and `ldms_grm_family` for the same
// reason, and the separation is load-bearing rather than tidy. An MPH relatedness matrix is `.grm.bin` plus
// `.grm.iid`, headed by an `int32` sample count and a `float32` sum of the SNP weights, over an unnormalised
// row-major upper triangle; a GCTA one is `.grm.bin`, `.grm.N.bin` and `.grm.id` over an already-normalised
// row-major lower triangle. Neither tool validates the other's layout and neither refuses it
// usefully: a GCTA
// bundle relabelled for MPH runs to exit 0 producing nothing, and an MPH bundle handed to GCTA is read without
// complaint -- `gcta --pca` on one returned eigenvalues of 11262.9, 7958.5 and 7701.4 against the true 2.78,
// 2.67 and 2.59, at exit 0. A shared backend name would let routing hand one estimator the other's bytes.
//
// `component_model` names the variance-component structure the *pipeline* plans for the estimator. It is what
// decides whether a row or a pair request may configure the LD- and MAF-stratified plan settings, in
// `resolveGctaMethodOptions`, `resolvePairRequests` and `resolveRelationships` — keyed on the planned model
// rather than on the matrix kind that carries it, so a second tool's stratified estimator can share one plan.
// It names the model the pipeline plans, not the number of components a native fit happens to report.
//
// `requires_complete_covariates` marks an estimator whose covariate interface neither reads nor reports a
// missing cell, so an incomplete file silently fits an arbitrary value. It gates the preparation rule.
def getMethodCapabilityContract() {
    return [
        required_fields: [
            'domain',
            'estimator_family',
            'input_backend',
            'component_model',
            'trait_support',
            'prevalence',
            'requires_complete_covariates',
            'citation_keys',
        ],
        queryable_fields: [
            'domain',
            'option_family',
            'matrix_kind',
            'endpoint_domain',
            'reference_family',
            'estimator_family',
            'input_backend',
            'component_model',
            'trait_support',
            'supports_covariates',
            'prevalence',
            'requires_complete_covariates',
            'citation_keys',
        ],
        estimator_families: [
            'whole_genome_regression',
            'mixed_linear_model',
            'reml',
            'moment_he',
            'pcgc',
            'ld_score_regression',
            'summary_tagging_regression',
            'meta_analysis',
            'meta_regression',
        ],
        input_backends: [
            'dense_grm',
            'ldms_grm_family',
            'ldak_kinship',
            'mph_grm',
            'mph_grm_family',
            'sparse_grm',
            'direct_plink_genotypes',
            'direct_plink1_genotypes',
            'summary_statistics',
        ],
        component_models: [
            'none',
            'single',
            'ld_maf_stratified',
            'tagging_bundle',
        ],
        genotype_bundles: [
            'plink',
            'plink1',
        ],
        prevalence_requirements: ['not_consumed', 'consumed', 'required'],
        trait_support_fields: ['quantitative', 'binary'],
    ]
}

// One registry owns selector domain, option family, matrix, prevalence, declared estimator capability and
// citation knowledge. Association entries also carry their GWASLab constructor mapping so the selector
// vocabulary cannot drift from it. Every capability is derived from what the wired module actually runs, not
// from the shape of the token, and `getMethodCapabilityContract()` fixes the vocabulary each one is written in.
//
// `prevalence` is the one typed scale contract. `consumed` accepts and uses a declared value, `required`
// additionally requires it at ingress, and `not_consumed` rejects it unless a downstream selected method uses
// it. HE and FASTHE convert binary estimates when population prevalence is declared; PCGC and FASTPCGC
// require it because their estimates are on the liability scale.
def getMethodRegistry() {
    return [
        fixed: [
            domain: 'summary_set',
            endpoint_domain: 'summary_statistics',
            estimator_family: 'meta_analysis',
            input_backend: 'summary_statistics',
            component_model: 'none',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'not_consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['gwaslab'],
        ],
        random: [
            domain: 'summary_set',
            endpoint_domain: 'summary_statistics',
            estimator_family: 'meta_analysis',
            input_backend: 'summary_statistics',
            component_model: 'none',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'not_consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['gwaslab'],
        ],
        re2: [
            domain: 'summary_set',
            endpoint_domain: 'summary_statistics',
            estimator_family: 'meta_analysis',
            input_backend: 'summary_statistics',
            component_model: 'none',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'not_consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['metasoft'],
        ],
        mrmega: [
            domain: 'summary_set',
            endpoint_domain: 'summary_statistics',
            estimator_family: 'meta_regression',
            input_backend: 'summary_statistics',
            component_model: 'none',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'not_consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['mrmega'],
        ],
        regenie: [
            domain: 'association',
            option_family: 'regenie',
            estimator_family: 'whole_genome_regression',
            input_backend: 'direct_plink_genotypes',
            component_model: 'none',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'not_consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['regenie'],
            mapping: [common: [
                snpid: 'ID',
                chrom: 'CHROM',
                pos: 'GENPOS',
                ea: 'ALLELE1',
                nea: 'ALLELE0',
                eaf: 'A1FREQ',
                n: 'N',
                beta: 'BETA',
                se: 'SE',
                mlog10p: 'LOG10P',
                readargs: [sep: ' '],
            ]],
        ],
        gcta_fastgwa: [
            domain: 'association',
            option_family: 'gcta',
            matrix_kind: 'gcta_sparse',
            estimator_family: 'mixed_linear_model',
            input_backend: 'sparse_grm',
            component_model: 'single',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'not_consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['gcta_fastgwa'],
            mapping: [common: [
                snpid: 'SNP',
                chrom: 'CHR',
                pos: 'POS',
                ea: 'A1',
                nea: 'A2',
                eaf: 'AF1',
                n: 'N',
                beta: 'BETA',
                se: 'SE',
                p: 'P',
            ]],
        ],
        ldak_kvik: [
            domain: 'association',
            option_family: 'ldak',
            estimator_family: 'mixed_linear_model',
            input_backend: 'direct_plink1_genotypes',
            component_model: 'single',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'not_consumed', sample: 'not_consumed'],
            requires_complete_covariates: true,
            citation_keys: ['ldak_kvik'],
            mapping: [
                common: [
                    snpid: 'Predictor',
                    chrom: 'Chromosome',
                    pos: 'Basepair',
                    ea: 'A1',
                    nea: 'A2',
                    eaf: 'EAF',
                    neff: 'N',
                    p: 'Wald_P',
                ],
                quantitative: [beta: 'Effect', se: 'SE'],
                binary: [beta: 'Approx_Log_OR', se: 'Approx_SE'],
            ],
        ],
        gcta_greml: [
            domain: 'heritability',
            option_family: 'gcta',
            matrix_kind: 'gcta_dense',
            estimator_family: 'reml',
            input_backend: 'dense_grm',
            component_model: 'single',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['gcta_greml'],
        ],
        gcta_greml_ldms: [
            domain: 'heritability',
            option_family: 'gcta',
            matrix_kind: 'gcta_ldms',
            estimator_family: 'reml',
            input_backend: 'ldms_grm_family',
            component_model: 'ld_maf_stratified',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['gcta_greml_ldms'],
        ],
        gcta_bivariate_reml: [
            domain: 'pairwise',
            endpoint_domain: 'analysis',
            option_family: 'gcta',
            matrix_kind: 'gcta_dense',
            estimator_family: 'reml',
            input_backend: 'dense_grm',
            component_model: 'single',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['gcta_bivariate_reml'],
        ],
        gcta_bivariate_reml_ldms: [
            domain: 'pairwise',
            endpoint_domain: 'analysis',
            option_family: 'gcta',
            matrix_kind: 'gcta_ldms',
            estimator_family: 'reml',
            input_backend: 'ldms_grm_family',
            component_model: 'ld_maf_stratified',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['gcta_bivariate_reml', 'gcta_greml_ldms'],
        ],
        gcta_bivariate_he: [
            domain: 'pairwise',
            endpoint_domain: 'analysis',
            option_family: 'gcta',
            matrix_kind: 'gcta_dense',
            estimator_family: 'moment_he',
            input_backend: 'dense_grm',
            component_model: 'single',
            trait_support: [quantitative: true, binary: false],
            supports_covariates: false,
            prevalence: [population: 'not_consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['gcta_hereg'],
        ],
        gcta_bivariate_he_ldms: [
            domain: 'pairwise',
            endpoint_domain: 'analysis',
            option_family: 'gcta',
            matrix_kind: 'gcta_ldms',
            estimator_family: 'moment_he',
            input_backend: 'ldms_grm_family',
            component_model: 'ld_maf_stratified',
            trait_support: [quantitative: true, binary: false],
            supports_covariates: false,
            prevalence: [population: 'not_consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['gcta_hereg', 'gcta_greml_ldms'],
        ],
        ldak_sumher: [
            domain: 'summary_unary',
            endpoint_domain: 'summary_statistics',
            option_family: 'ldak',
            reference_family: 'ldak',
            estimator_family: 'summary_tagging_regression',
            input_backend: 'summary_statistics',
            component_model: 'tagging_bundle',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'consumed', sample: 'consumed'],
            requires_complete_covariates: false,
            citation_keys: ['ldak_sumstats'],
        ],
        ldak_sumcors: [
            domain: 'pairwise',
            endpoint_domain: 'summary_statistics',
            option_family: 'ldak',
            reference_family: 'ldak',
            estimator_family: 'summary_tagging_regression',
            input_backend: 'summary_statistics',
            component_model: 'tagging_bundle',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'consumed', sample: 'consumed'],
            requires_complete_covariates: false,
            citation_keys: ['ldak_sumstats'],
        ],
        ldsc_h2: [
            domain: 'summary_unary',
            endpoint_domain: 'summary_statistics',
            option_family: 'ldsc',
            reference_family: 'ldsc',
            estimator_family: 'ld_score_regression',
            input_backend: 'summary_statistics',
            component_model: 'single',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'consumed', sample: 'consumed'],
            requires_complete_covariates: false,
            citation_keys: ['ldsc'],
        ],
        ldsc_rg: [
            domain: 'pairwise',
            endpoint_domain: 'summary_statistics',
            option_family: 'ldsc',
            reference_family: 'ldsc',
            estimator_family: 'ld_score_regression',
            input_backend: 'summary_statistics',
            component_model: 'single',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'consumed', sample: 'consumed'],
            requires_complete_covariates: false,
            citation_keys: ['ldsc'],
        ],
        ldak_reml: [
            domain: 'heritability',
            option_family: 'ldak',
            matrix_kind: 'ldak_kinship',
            estimator_family: 'reml',
            input_backend: 'ldak_kinship',
            component_model: 'single',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'consumed', sample: 'not_consumed'],
            requires_complete_covariates: true,
            citation_keys: ['ldak'],
        ],
        ldak_he: [
            domain: 'heritability',
            option_family: 'ldak',
            matrix_kind: 'ldak_kinship',
            estimator_family: 'moment_he',
            input_backend: 'ldak_kinship',
            component_model: 'single',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'consumed', sample: 'not_consumed'],
            requires_complete_covariates: true,
            citation_keys: ['ldak'],
        ],
        ldak_pcgc: [
            domain: 'heritability',
            option_family: 'ldak',
            matrix_kind: 'ldak_kinship',
            estimator_family: 'pcgc',
            input_backend: 'ldak_kinship',
            component_model: 'single',
            trait_support: [quantitative: false, binary: true],
            prevalence: [population: 'required', sample: 'not_consumed'],
            requires_complete_covariates: true,
            citation_keys: ['ldak'],
        ],
        ldak_fast_he: [
            domain: 'heritability',
            option_family: 'ldak',
            estimator_family: 'moment_he',
            input_backend: 'direct_plink1_genotypes',
            component_model: 'single',
            trait_support: [quantitative: true, binary: true],
            prevalence: [population: 'consumed', sample: 'not_consumed'],
            requires_complete_covariates: true,
            citation_keys: ['ldak', 'rhe_mc'],
        ],
        ldak_fast_pcgc: [
            domain: 'heritability',
            option_family: 'ldak',
            estimator_family: 'pcgc',
            input_backend: 'direct_plink1_genotypes',
            component_model: 'single',
            trait_support: [quantitative: false, binary: true],
            prevalence: [population: 'required', sample: 'not_consumed'],
            requires_complete_covariates: true,
            citation_keys: ['ldak', 'rhe_mc'],
        ],
        mph_reml: [
            domain: 'heritability',
            option_family: 'mph',
            matrix_kind: 'mph_dense',
            estimator_family: 'reml',
            input_backend: 'mph_grm',
            component_model: 'single',
            trait_support: [quantitative: true, binary: false],
            prevalence: [population: 'not_consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['mph'],
        ],
        mph_reml_ldms: [
            domain: 'heritability',
            option_family: 'mph',
            matrix_kind: 'mph_ldms',
            estimator_family: 'reml',
            input_backend: 'mph_grm_family',
            component_model: 'ld_maf_stratified',
            trait_support: [quantitative: true, binary: false],
            prevalence: [population: 'not_consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['mph', 'gcta_greml_ldms'],
        ],
        mph_bivariate_reml: [
            domain: 'pairwise',
            endpoint_domain: 'analysis',
            option_family: 'mph',
            matrix_kind: 'mph_dense',
            estimator_family: 'reml',
            input_backend: 'mph_grm',
            component_model: 'single',
            trait_support: [quantitative: true, binary: false],
            prevalence: [population: 'not_consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['mph'],
        ],
        mph_bivariate_reml_ldms: [
            domain: 'pairwise',
            endpoint_domain: 'analysis',
            option_family: 'mph',
            matrix_kind: 'mph_ldms',
            estimator_family: 'reml',
            input_backend: 'mph_grm_family',
            component_model: 'ld_maf_stratified',
            trait_support: [quantitative: true, binary: false],
            prevalence: [population: 'not_consumed', sample: 'not_consumed'],
            requires_complete_covariates: false,
            citation_keys: ['mph', 'gcta_greml_ldms'],
        ],
    ]
}

// Preserve the established public mapping shape while projecting it from the unified registry.
def getAssociationColumnMappings() {
    return getMethodRegistry()
        .findAll { _token, entry -> entry.domain == 'association' }
        .collectEntries { token, entry -> [(token): entry.mapping] }
}

def getAssociationColumnMappingJson(method, is_binary) {
    def entry = getAssociationColumnMappings()[method]
    if (!entry) {
        error("[nf-core/gwas] ERROR: no GWASLab column mapping is registered for association method '${method}'")
    }
    return groovy.json.JsonOutput.toJson(entry.common + (entry[is_binary ? 'binary' : 'quantitative'] ?: [:]))
}

// Project the stable capability interface without exposing the association mapping implementation.
def getMethodCapabilities() {
    return getMethodRegistry().collectEntries { token, entry ->
        [(token): entry.findAll { name, _value -> name != 'mapping' }]
    }
}

// Ask the registry which methods hold a capability. Every caller that used to carry its own method-name list
// goes through here, so adding an entry extends the answer instead of leaving one list quietly behind.
def getMethodTokensWithCapabilities(required) {
    def contract = getMethodCapabilityContract()
    def unknown = required.keySet().findAll { field -> !(field in contract.queryable_fields) }
    if (unknown) {
        error("[nf-core/gwas] ERROR: method capability selection uses unregistered field '${unknown.first()}'")
    }
    return getMethodCapabilities()
        .findAll { _token, details -> required.every { field, value -> details[field] == value } }
        .keySet()
        .toList()
}

def getMethodCapability(method, field) {
    def details = getMethodCapabilities()[method]
    if (!details) {
        error("[nf-core/gwas] ERROR: no capability is registered for method '${method}'")
    }
    if (!details.containsKey(field)) {
        error("[nf-core/gwas] ERROR: method '${method}' declares no capability '${field}'")
    }
    return details[field]
}

// A direct-genotype backend names the genotype representation the estimator's own executable reads.
// `direct_plink_genotypes` is format-polymorphic: REGENIE picks `--bfile` or `--pfile` from the primary
// staged extension. `direct_plink1_genotypes` is BED/BIM/FAM only, which is what every LDAK executable reads.
def getDirectGenotypeBackendContract() {
    return [
        direct_plink_genotypes: [genotype_bundle: 'plink'],
        direct_plink1_genotypes: [genotype_bundle: 'plink1'],
    ]
}

// The closed set of reusable relatedness intermediates, each declaring the backend its consumers read and the
// genotype representation its builder needs. Iteration order is the order matrix kinds are requested in, so it
// is part of the contract rather than an incidental map layout.
def getMatrixKindContract() {
    return [
        gcta_dense: [input_backend: 'dense_grm', genotype_bundle: 'plink'],
        gcta_ldms: [input_backend: 'ldms_grm_family', genotype_bundle: 'plink1'],
        gcta_sparse: [input_backend: 'sparse_grm', genotype_bundle: 'plink'],
        ldak_kinship: [input_backend: 'ldak_kinship', genotype_bundle: 'plink1'],
        mph_dense: [input_backend: 'mph_grm', genotype_bundle: 'plink1'],
        mph_ldms: [input_backend: 'mph_grm_family', genotype_bundle: 'plink1'],
    ]
}

// `null` for a summary estimator, which consumes no genotypes at all; otherwise the representation the
// estimator itself or its matrix builder reads.
def getMethodGenotypeBundle(method) {
    def details = getMethodCapabilities()[method]
    if (!details) {
        error("[nf-core/gwas] ERROR: no capability is registered for method '${method}'")
    }
    if (details.input_backend == 'summary_statistics') {
        return null
    }
    // Every caller of this function iterates the whole registry, so an undescribed backend must fail by name
    // here rather than as a null dereference that reports nothing about which entry is at fault.
    def contract = details.matrix_kind
        ? getMatrixKindContract()[details.matrix_kind]
        : getDirectGenotypeBackendContract()[details.input_backend]
    if (!contract) {
        error("[nf-core/gwas] ERROR: method '${method}' declares ${details.matrix_kind ? "matrix kind '${details.matrix_kind}'" : "input backend '${details.input_backend}'"}, which no genotype-bundle contract describes")
    }
    return contract.genotype_bundle
}

// The PLINK 1 compatibility bundle is a genotype-representation requirement, not a scientific one: a method
// needs the derivative when its direct backend or its matrix kind declares the `plink1` bundle. Deriving it
// from the declared bundle rather than from `option_family == 'ldak'` is what lets a non-LDAK direct-genotype
// estimator join the answer without editing this predicate.
def getPlink1GenotypeMethodTokens() {
    return getMethodCapabilities()
        .keySet()
        .findAll { method -> getMethodGenotypeBundle(method) == 'plink1' }
        .toList()
}

def getAssociationMethodTokens() {
    return getAssociationColumnMappings().keySet().toList()
}

def getHeritabilityMethodTokens() {
    return getMethodCapabilities()
        .findAll { _token, details -> details.domain == 'heritability' }
        .keySet()
        .toList()
}

def getSummaryUnaryMethodTokens() {
    return getMethodCapabilities()
        .findAll { _token, details -> details.domain == 'summary_unary' }
        .keySet()
        .toList()
}

def getRelationshipMethodTokens() {
    return getMethodCapabilities()
        .findAll { _token, details -> details.domain == 'pairwise' }
        .keySet()
        .toList()
}

def getSummarySetMethodTokens() {
    return getMethodTokensWithCapabilities([domain: 'summary_set'])
}
