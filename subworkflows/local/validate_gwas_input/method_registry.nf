// The closed vocabularies every registry entry is written against. Routing and validation ask the registry
// which methods hold a capability instead of repeating a method-name list that silently goes stale when a
// method is added, so a new entry that omits or misspells a capability must fail the registry contract test
// rather than quietly fall out of a hardcoded list.
//
// `input_backend` names the representation the estimator genuinely consumes. A direct-genotype estimator is a
// distinct backend and must never be given a fabricated matrix kind to make it look like a GRM route. The
// direct backends are split by genotype representation, not by tool: `direct_plink1_genotypes` names an
// executable that reads BED/BIM/FAM only, `direct_plink_genotypes` one that selects the native flag from the
// staged primary extension. Two tools' relatedness matrices never share a backend name even when their file
// extensions coincide, because a shared name would let routing hand one estimator the other's bytes.
//
// `component_model` names the variance-component structure the *pipeline* plans for the estimator. It is what
// decides whether a row may configure LD- and MAF-stratified plan settings, keyed on the planned model rather
// than on the matrix kind that carries it, so a second tool's stratified estimator can share one plan.
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
// citation knowledge, so the selector vocabulary cannot drift from the wired routes. Every capability is
// derived from what the wired module actually runs, not from the shape of the token, and
// `getMethodCapabilityContract()` fixes the vocabulary each one is written in.
//
// `prevalence` is the one typed scale contract. `consumed` accepts and uses a declared value, `required`
// additionally requires it at ingress, and `not_consumed` rejects it unless a downstream selected method uses
// it.
//
// No analysis method is wired yet: the registry is empty, so both selector vocabularies are empty and every
// populated selector cell is rejected at ingress. Each analysis route adds its entry here when it is wired.
def getMethodRegistry() {
    return [:]
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
// `direct_plink_genotypes` is format-polymorphic: the estimator picks its PLINK 1 or PLINK 2 input flag from the
// primary staged extension. `direct_plink1_genotypes` is BED/BIM/FAM only.
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
    // No matrix-backed method is registered yet, so no reusable relatedness intermediate exists.
    return [:]
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
// from the declared bundle rather than from the option family is what lets any direct-genotype estimator join
// the answer without editing this predicate.
def getPlink1GenotypeMethodTokens() {
    return getMethodCapabilities()
        .keySet()
        .findAll { method -> getMethodGenotypeBundle(method) == 'plink1' }
        .toList()
}

def getAssociationMethodTokens() {
    return getMethodRegistry()
        .findAll { _token, entry -> entry.domain == 'association' }
        .keySet()
        .toList()
}

def getHeritabilityMethodTokens() {
    return getMethodCapabilities()
        .findAll { _token, details -> details.domain == 'heritability' }
        .keySet()
        .toList()
}
