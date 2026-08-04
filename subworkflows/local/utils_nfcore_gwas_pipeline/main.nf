// Pipeline-specific initialisation, relational input validation and completion utilities.
// These workflows contain no version-producing processes, so they emit no versions.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// SUBWORKFLOW: Consisting entirely of nf-core/modules
include { UTILS_NFSCHEMA_PLUGIN   } from '../../nf-core/utils_nfschema_plugin'
include { completionEmail         } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary       } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'

// PLUGIN
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { samplesheetToList       } from 'plugin/nf-schema'
include { paramsHelp              } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version // channel: val(version)
    validate_params // channel: val(validate_params)
    monochrome_logs // channel: val(monochrome_logs)
    nextflow_cli_args // channel: val(nextflow_cli_args)
    outdir // channel: val(outdir)
    cohort_manifest // channel: val(cohort_manifest)
    analysis_manifest // channel: val(analysis_manifest)
    method_options // channel: val(method_options)
    help // channel: val(help)
    help_full // channel: val(help_full)
    show_hidden // channel: val(show_hidden)

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/gwas ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/', '')}" }.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/gwas/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --cohort_manifest cohorts.csv --analysis_manifest analyses.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Validate both manifests as complete lists before channel construction, so cross-row and
    // cross-manifest failures cannot submit downstream tasks.
    //
    if (!cohort_manifest) {
        error("[nf-core/gwas] ERROR: --analysis_manifest '${analysis_manifest ?: ''}' was supplied but --cohort_manifest is missing; linked manifests require both")
    }
    if (!analysis_manifest) {
        error("[nf-core/gwas] ERROR: --cohort_manifest '${cohort_manifest}' was supplied but --analysis_manifest is missing; linked manifests require both")
    }
    def cohort_schema = "${projectDir}/assets/schema_cohort_manifest.json"
    def analysis_schema = "${projectDir}/assets/schema_analysis_manifest.json"
    validateSamplesheetHeader(cohort_manifest, cohort_schema, 'Cohort manifest')
    validateSamplesheetHeader(analysis_manifest, analysis_schema, 'Analysis manifest')
    def ch_analyses = channel.fromList(
        validateRelationalInput(
            samplesheetToList(cohort_manifest, cohort_schema),
            samplesheetToList(analysis_manifest, analysis_schema),
            cohort_manifest,
            analysis_manifest,
            cohort_schema,
            analysis_schema,
            method_options,
        )
    )

    emit:
    analyses = ch_analyses // channel: [ val(meta), [ path(genotype_file), ... ], path(phenotype), path(quant_covariates), path(cat_covariates), path(kvik_extract), path(ldak_weights) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email // channel: val(email)
    email_on_fail // channel: val(email_on_fail)
    plaintext_email // channel: val(plaintext_email)
    outdir // channel: val(outdir)
    monochrome_logs // channel: val(monochrome_logs)
    multiqc_report // channel: [ [ path(report) ] ]

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// The GWASLab constructor column mapping for each association method.
//
// Explicit mappings rather than GWASLab format names, for every method and not only the ones with no
// format-book entry: LDAK-KVIK has no entry at all, the pinned GWASLab version differs from the version
// whose format book was inspected, and PLINK 2's own entry declares a `#` comment character, which would
// make pandas read the `#CHROM` header line as a comment. Four short maps cost less than four unverified
// assumptions.
//
// The keys are `gwaslab.Sumstats` constructor arguments and the values are the source column headers the
// programme actually writes; the component splats the serialised map straight into that constructor. A key
// that is not a constructor argument is forwarded to `pandas.read_table` instead, which is how `readargs`
// carries a non-tab separator.
//
// `common` applies to every analysis. `quantitative` and `binary` are merged over it where the programme
// renames its effect columns by trait type; a method that does not need them omits both.
//
// Registering a new association route is exactly this: add an entry. The accepted `association_methods`
// vocabulary is derived from these keys, so a route cannot be selectable without a mapping and a mapping
// cannot be orphaned.
//
def getAssociationColumnMappings() {
    return [
        plink2: [common: [
            snpid: 'ID',
            chrom: '#CHROM',
            pos: 'POS',
            ea: 'A1',
            nea: 'REF',
            eaf: 'A1_FREQ',
            n: 'OBS_CT',
            beta: 'BETA',
            se: 'SE',
            p: 'P',
        ]],
        regenie: [common: [
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
        gcta_fastgwa: [common: [
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
        ldak_kvik: [
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
    ]
}

//
// The column mapping for one analysis, serialised as the JSON object the component takes as its
// `input_format`. Serialised here rather than written out by hand so that quoting is the JSON library's
// problem rather than a reviewer's.
//
def getAssociationColumnMappingJson(method, is_binary) {
    def entry = getAssociationColumnMappings()[method]
    if (!entry) {
        error("[nf-core/gwas] ERROR: no GWASLab column mapping is registered for association method '${method}'")
    }
    return groovy.json.JsonOutput.toJson(entry.common + (entry[is_binary ? 'binary' : 'quantitative'] ?: [:]))
}

//
// The optional GWASLab reference resources, keyed on genome build and on nothing else. The samplesheet's
// `genome_build` column is the only real heterogeneity in this set — an rsID VCF is a property of the
// coordinate system, not of a cohort or an ancestry — so the reference selection reduces to this.
//
// A build with no configured resource yields empty lists, which stage nothing and reach the component as
// an absent reference; harmonisation then standardises without one rather than failing. Index sidecars are
// read by convention rather than exposed as six further parameters, and a missing one is a hard error
// naming the file, because a silently absent index would fail later inside GWASLab.
//
def getGwaslabReferences() {
    def sidecar = { resource, suffixes ->
        if (!resource) {
            return []
        }
        def resolved = suffixes.collect { suffix -> file("${resource}${suffix}") }.find { candidate -> candidate.exists() }
        if (!resolved) {
            error("[nf-core/gwas] ERROR: no index found for GWASLab reference '${resource}', expected one of ${suffixes.collect { suffix -> "'${resource}${suffix}'" }.join(', ')}")
        }
        return resolved
    }
    def resources = { fasta, rsid_vcf, strand_vcf ->
        [
            fasta: fasta ? file(fasta, checkIfExists: true) : [],
            fasta_index: sidecar.call(fasta, ['.fai']),
            rsid_vcf: rsid_vcf ? file(rsid_vcf, checkIfExists: true) : [],
            rsid_vcf_index: sidecar.call(rsid_vcf, ['.tbi', '.csi']),
            strand_vcf: strand_vcf ? file(strand_vcf, checkIfExists: true) : [],
            strand_vcf_index: sidecar.call(strand_vcf, ['.tbi', '.csi']),
        ]
    }
    return [
        GRCh37: resources.call(params.gwaslab_reference_fasta_grch37, params.gwaslab_rsid_vcf_grch37, params.gwaslab_strand_vcf_grch37),
        GRCh38: resources.call(params.gwaslab_reference_fasta_grch38, params.gwaslab_rsid_vcf_grch38, params.gwaslab_strand_vcf_grch38),
    ]
}

//
// Accepted association method tokens for the `association_methods` column. Derived from the column-mapping
// registry: every association result is harmonised, so a route that is selectable without a registered
// mapping is a route that fails at harmonisation time. The literal order of that map is preserved, so the
// vocabulary this reports in a validation error is unchanged.
//
def getAssociationMethodTokens() {
    return getAssociationColumnMappings().keySet().toList()
}

//
// Accepted heritability method tokens for the `heritability_methods` column. The three LDAK
// estimators are separate tokens rather than one `ldak` token behind a sub-selector, so a single
// analysis can request all three.
//
def getHeritabilityMethodTokens() {
    return ['gcta_greml', 'gcta_greml_ldms', 'ldak_reml', 'ldak_he', 'ldak_pcgc']
}

//
// One declared value rendered as text, so that two values of the same meaning render alike.
//
// A declared value is a quantity a researcher stated in structured method options. JSON numbers may
// reach Groovy as an Integer or BigDecimal while equivalent values can be represented as strings.
// Comparing them unrendered would give one matrix two keys and build it twice.
//
// Numbers, and strings that spell a number, both render through BigDecimal with trailing zeros stripped, so
// 4, 4.0, '4' and '4.00' are one value and -0.25 and '-.25' are another. `BigDecimal` is constructed from
// the value's text rather than from a double, which would render 0.1 as 0.1000000000000000055511151231257827
// and split a CSV `0.1` from a schema-default `0.1`. Nulls render distinctly from the empty string, because
// "the researcher declared nothing" should not be silently conflated with a literal zero.
//
// Identifiers are rendered by canonicaliseIdentifier below instead, and never by this function: numeric
// equivalence is meaningless between names and collapsing them would over-share a matrix.
//
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

//
// One identifier rendered as text, verbatim.
//
// An identifier is a name the researcher chose or the pipeline assigned — a cohort id, a genotype format, a
// genotype file name, a matrix kind — and two names are the same only when they are spelled the same. Passing
// one through canonicaliseDeclaredValue instead would read numeric-looking names as numbers and collapse the
// cohorts `01` and `1` (likewise `7`, `7.0`, `+7`, and `1e3` with `1000`) onto one reuse key, so the analyses
// of one cohort would silently be estimated against the other cohort's relatedness matrix. Every other seam
// that decides cohort identity — the preflight's genotype agreement check, PREPARE_COHORT_GENOTYPES'
// deduplication — compares these names textually, and the key has to agree with them or the two layers
// disagree about how many cohorts a run has.
//
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

//
// The reuse key of one relatedness matrix: a digest over the cohort identity and every declared input that
// changes the matrix, and over nothing else.
//
// Identity and settings are taken as two arguments rather than as one map because the two are rendered by
// different rules — identifiers verbatim, declared values numerically — and a single map would leave the
// choice to be made per component name, which is exactly the decision a route registering a new setting
// should not have to make.
//
// Content-derived rather than arrival-ordered, so that adding a row to the samplesheet does not renumber the
// matrices that were already built and invalidate resume. Truncated to twelve hex characters because it is a
// path segment a researcher has to read and quote; twelve characters is 48 bits, which is far beyond
// collision range for the tens of matrices one run can declare.
//
def buildRelatednessMatrixKey(identity, settings) {
    def rendered = identity.collectEntries { name, value -> [(name): canonicaliseIdentifier(value)] } + [settings: canonicaliseDeclaredValue(settings)]
    def canonical = rendered
        .sort { entry -> entry.key }
        .collect { name, text -> "${name}=${text}" }
        .join('\n')
    return java.security.MessageDigest
        .getInstance('SHA-256')
        .digest(canonical.getBytes('UTF-8'))
        .encodeHex()
        .toString()
        .substring(0, 12)
}

//
// The relatedness matrices one analysis unit needs, by the methods it selected.
//
// A matrix kind names the tool and the matrix family together, because cross-tool reuse is not possible —
// GCTA and LDAK matrices are different formats produced by different algorithms — so the kind is itself a key
// component and two tools can never collide on one key.
//
// Registering a route is exactly this: add its method token here and its construction settings in
// getRelatednessMatrixSettings below. Adding a route's settings anywhere else silently makes two different
// matrices share one key.
//
def getRelatednessMatrixKinds(meta) {
    def kinds = []
    if ('gcta_greml' in meta.heritability_methods) {
        kinds << 'gcta_dense'
    }
    if ('gcta_greml_ldms' in meta.heritability_methods) {
        kinds << 'gcta_ldms'
    }
    if ('gcta_fastgwa' in meta.association_methods) {
        kinds << 'gcta_sparse'
    }
    if (meta.heritability_methods.any { method -> method in ['ldak_reml', 'ldak_he', 'ldak_pcgc'] }) {
        kinds << 'ldak_kinship'
    }
    return kinds.unique()
}

//
// The portable identity of an optional LDAK weights file.
//
// Only the file's bytes decide whether two analyses can reuse one matrix. A path or basename would make
// published keys machine-specific, would fail to reuse identical files copied elsewhere and would collapse
// different files sharing one name. The actual Path therefore stays on the build tuple outside this map.
//
def getLdakWeightsIdentity(weights_file, weights_policy = 'equal') {
    if (!weights_file) {
        return [mode: weights_policy]
    }
    def weights_path = weights_file instanceof java.nio.file.Path ? weights_file : weights_file.toPath()
    def digest = java.security.MessageDigest.getInstance('SHA-256')
    java.nio.file.Files
        .newInputStream(weights_path)
        .withCloseable { input ->
            input.eachByte(8192) { buffer, count ->
                digest.update(buffer, 0, count)
            }
        }
    return [mode: 'provided', sha256: digest.digest().encodeHex().toString()]
}

//
// The declared construction settings of one matrix kind: the inputs that change the matrix itself. Every
// GCTA matrix family reads its validated per-analysis family map. Dense matrices key the MAF threshold and
// the content identity of an optional extraction resource; LDMS keys all three stratification settings;
// sparse matrices key their cutoff. Estimator-only options are deliberately absent. LDAK's model, power,
// weights and relatedness filter are key components: filtered and all-sample matrices are scientifically
// different artefacts and therefore cannot share a deduplicated matrix identity.
//
// An `if` chain rather than a `switch`: `nextflow lint` aborts on any `switch` statement it is given
// (`ERROR ~ begin N, end N+1, length N`), so the construct cannot appear in this repository at all.
//
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

//
// One analysis unit's request for one matrix: its reuse key, and the construction inputs that are
// deliberately outside the key and so have to be reconciled across the rows that share it.
//
// The genotypes enter the key as sorted basenames rather than as paths: absolute paths differ per machine
// and the digest is a published path segment, so a path-derived digest would make every published tree
// machine-specific.
//
// `parts` is deliberately absent because GCTA's --make-grm-part count is run/profile configuration. It
// changes task partitioning but not matrix content and therefore belongs neither in analysis metadata nor
// in the reuse key.
//
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


//
// The three mutually exclusive genotype groups, each mapped to the columns that make it complete.
//
def getGenotypeGroups() {
    return [
        plink2: ['pgen', 'psam', 'pvar'],
        plink1: ['bed', 'bim', 'fam'],
        vcf: ['vcf'],
    ]
}


//
// The samplesheet columns nf-schema emits positionally after the meta map, in the order it emits
// them. Derived from the schema rather than restated here: a column added there shifts every later
// slot, and a hand-maintained copy would silently start reading cells from the wrong column. A
// column is positional exactly when it declares no `meta` key, which today makes the list identical
// to the columns carrying files.
//
def getSamplesheetPositionalColumns(schema) {
    def properties = new groovy.json.JsonSlurper().parseText(file(schema).text).items.properties
    return properties.findAll { _column, definition -> !definition.containsKey('meta') }.keySet().toList()
}

//
// Each CSV contract requires its complete header, including optional cells. JSON Schema deliberately
// does not require optional cells and nf-schema would therefore accept a CSV that omitted their
// columns entirely. Compare the actual header before nf-schema injects defaults. Column order is not
// significant, but missing, unexpected or repeated names are all contract errors.
//
def validateSamplesheetHeader(samplesheet, schema, role) {
    def expected = new groovy.json.JsonSlurper().parseText(file(schema).text).items.properties.keySet().toList()
    def header_line = file(samplesheet).readLines().find { line -> line.trim() }
    def observed = header_line
        ? header_line.split(',', -1).collect { column -> column.trim().replaceAll(/^"|"$/, '') }
        : []
    def missing = expected.findAll { column -> !observed.contains(column) }
    def unexpected = observed.findAll { column -> !expected.contains(column) }.unique()
    def repeated = observed.countBy { column -> column }.findAll { _column, count -> count > 1 }.keySet().toList()

    def problems = []
    if (missing) {
        problems << "missing column headers: ${missing.collect { column -> "'${column}'" }.join(', ')}"
    }
    if (unexpected) {
        problems << "unexpected column headers: ${unexpected.collect { column -> "'${column}'" }.join(', ')}"
    }
    if (repeated) {
        problems << "repeated column headers: ${repeated.collect { column -> "'${column}'" }.join(', ')}"
    }
    if (problems) {
        error("[nf-core/gwas] ERROR: ${role} '${samplesheet}' header row 1 does not match the mandatory ${expected.size()}-column input contract.\n\n  - ${problems.join('\n  - ')}\n")
    }
}

//
// An empty samplesheet cell never reaches the meta map as null: nf-schema drops the key and its
// converter substitutes an empty list in its place. A populated cell, meanwhile, may arrive as a
// number whose Groovy truth is false — `0` is both a legal PLINK control code and a legal MAF bin
// edge. So "was this column populated?" is asked through this, and its answer compared against
// null rather than taken as a truth value.
//
def normaliseCellValue(value) {
    if (value == null || (value instanceof Collection && value.isEmpty())) {
        return null
    }
    return value.toString().trim() ? value : null
}

//
// Split a comma-delimited method selector cell into its tokens.
//
def tokenizeMethodSelector(selector) {
    return selector ? selector.toString().tokenize(',').collect { token -> token.trim() }.findAll { token -> token } : []
}

//
// The method routes a row selects, as the membership tests the conditional columns are gated on.
// Each is membership of a named token set rather than "the selector list is non-empty", so a row
// naming only an unrecognised method is still treated as selecting nothing.
//
def getMethodRoutes(association_methods, heritability_methods) {
    return [
        runs_heritability: heritability_methods.any { method -> method in getHeritabilityMethodTokens() },
        consumes_population_prevalence: heritability_methods.any { method ->
            method in ['gcta_greml', 'gcta_greml_ldms', 'ldak_reml', 'ldak_pcgc']
        },
        runs_ldak_kvik: association_methods.contains('ldak_kvik'),
        runs_ldak_heritability: heritability_methods.any { method -> method in ['ldak_reml', 'ldak_he', 'ldak_pcgc'] },
        runs_ldak_pcgc: heritability_methods.contains('ldak_pcgc'),
        runs_gcta: association_methods.contains('gcta_fastgwa') || heritability_methods.any { method -> method in ['gcta_greml', 'gcta_greml_ldms'] },
        runs_gcta_fastgwa: association_methods.contains('gcta_fastgwa'),
        runs_greml_ldms: heritability_methods.contains('gcta_greml_ldms'),
    ]
}

//
// Analysis-manifest trait semantics resolved once before validation and canonical tuple construction.
//
def getAnalysisSettings(meta) {
    return [
        population_prevalence: normaliseCellValue(meta.population_prevalence),
        case_value: normaliseCellValue(meta.case_value),
        control_value: normaliseCellValue(meta.control_value),
    ]
}

//
// Method selectors: every token must belong to its column's vocabulary and may appear only once.
// These are resolved before anything else because membership of them is what makes nine further
// columns meaningful or meaningless.
//
def validateMethodSelectors(association_methods, heritability_methods, reject) {
    [
        [column: 'association_methods', methods: association_methods, vocabulary: getAssociationMethodTokens()],
        [column: 'heritability_methods', methods: heritability_methods, vocabulary: getHeritabilityMethodTokens()],
    ].each { selector ->
        def unknown = selector.methods.findAll { method -> !selector.vocabulary.contains(method) }.unique()
        if (unknown) {
            reject.call(selector.column, "unknown method${unknown.size() > 1 ? 's' : ''} ${unknown.collect { method -> "'${method}'" }.join(', ')}, accepted values are ${selector.vocabulary.collect { method -> "'${method}'" }.join(', ')}")
        }
        def repeated = selector.methods.countBy { method -> method }.findAll { _method, count -> count > 1 }.keySet()
        if (repeated) {
            reject.call(selector.column, "method${repeated.size() > 1 ? 's' : ''} ${repeated.collect { method -> "'${method}'" }.join(', ')} listed more than once")
        }
    }
    // Neither selector alone is required, so an empty row is only detectable from both together and
    // is attributed to both.
    if (!association_methods && !heritability_methods) {
        reject.call(['association_methods', 'heritability_methods'], "row selects no method, populate one of them or remove the row")
    }
}

//
// Genotype groups: exactly one, populated completely. Returns the name of the one complete group,
// or null when the row does not have one.
//
def validateGenotypeGroup(cells, reject) {
    def groups = getGenotypeGroups()
    def populated_groups = groups.findAll { _name, columns -> columns.any { column -> cells[column] } }
    if (!populated_groups) {
        // Nothing on the row points at a genotype column, so every genotype column is equally
        // implicated and naming one of them would be arbitrary.
        reject.call(groups.values().flatten(), "no genotype group is populated, supply exactly one of pgen/psam/pvar, bed/bim/fam or vcf")
    }
    else if (populated_groups.size() > 1) {
        populated_groups
            .keySet()
            .toList()
            .tail()
            .each { name ->
                def populated_column = groups[name].find { column -> cells[column] }
                reject.call(populated_column, "a second genotype group is populated on this row, supply exactly one of pgen/psam/pvar, bed/bim/fam or vcf")
            }
    }
    else {
        def group_name = populated_groups.keySet().first()
        groups[group_name]
            .findAll { column -> !cells[column] }
            .each { column ->
                reject.call(column, "genotype group '${group_name}' is only partly populated, all of ${groups[group_name].join(', ')} are required together")
            }
    }
    return populated_groups.size() == 1 ? populated_groups.keySet().first() : null
}


//
// Trait type drives the case, control and prevalence columns. It is the canonical trait
// classification and is never inferred from the values or from prevalence.
//
def validateTraitColumns(is_binary, settings, reject) {
    if (is_binary) {
        if (settings.case_value == null) {
            reject.call('case_value', "a binary trait must declare the value used for cases in the phenotype file")
        }
        if (settings.control_value == null) {
            reject.call('control_value', "a binary trait must declare the value used for controls in the phenotype file")
        }
        if (settings.case_value != null && settings.control_value != null && settings.case_value.toString() == settings.control_value.toString()) {
            reject.call(['case_value', 'control_value'], "binary case_value and control_value must be distinct source codes")
        }
    }
    else {
        if (settings.case_value != null) {
            reject.call('case_value', "case_value has no meaning on a quantitative trait, remove it or set trait_type to 'binary'")
        }
        if (settings.control_value != null) {
            reject.call('control_value', "control_value has no meaning on a quantitative trait, remove it or set trait_type to 'binary'")
        }
        if (settings.population_prevalence != null) {
            reject.call('population_prevalence', "prevalence has no meaning on a quantitative trait, remove it or set trait_type to 'binary'")
        }
    }
}

//
// Population prevalence belongs to heritability analyses and is mandatory for PCGC.
//
def validateMethodConditionedColumns(settings, routes, reject) {
    if (settings.population_prevalence != null && !routes.consumes_population_prevalence) {
        reject.call('population_prevalence', "none of the selected estimators consumes it; select GCTA GREML, GCTA GREML-LDMS, LDAK REML or LDAK PCGC, or remove the prevalence")
    }
    if (settings.population_prevalence == null && routes.runs_ldak_pcgc) {
        reject.call('population_prevalence', "'ldak_pcgc' always estimates on the liability scale and requires a population prevalence")
    }
}


//
// Portable content identity for a stageable method resource. Paths are deliberately excluded: copied
// resources with identical bytes reuse one matrix, while changed bytes under one basename do not collide.
//
def getMethodResourceIdentity(resource) {
    if (!resource) {
        return [mode: 'all']
    }
    def resource_path = resource instanceof java.nio.file.Path ? resource : resource.toPath()
    def digest = java.security.MessageDigest.getInstance('SHA-256')
    java.nio.file.Files
        .newInputStream(resource_path)
        .withCloseable { input ->
            input.eachByte(8192) { buffer, count ->
                digest.update(buffer, 0, count)
            }
        }
    return [mode: 'file', sha256: digest.digest().encodeHex().toString()]
}

//
// Parse and validate the advanced method-options document before any channels are constructed. The root is
// keyed by analysis_id and each value is a method-family map. All GCTA matrix construction and estimator
// settings share this typed family vocabulary so scientific inputs cannot bypass validation or reuse
// identity.
//
def validateMethodOptions(method_options, analysis_rows) {
    def defaults = [
        gcta: [
            grm_maf: null,
            grm_extract: [],
            reml_no_constrain: false,
            sparse_cutoff: 0.05,
            ld_score_region_kb: 200,
            ld_bins: 4,
            ldms_maf_edges: [0, 0.01, 0.05, 0.2, 0.5],
        ],
        ldak: [
            model: 'human_default',
            power: -0.25,
            weights_policy: 'equal',
            weights: [],
            relatedness_filter: false,
            kvik_step1_subset: 'all',
            predictor_extract: [],
        ],
        regenie: [
            step1_bsize: 1000,
            firth: true,
            firth_approx: true,
            firth_p_threshold: 0.01,
            min_mac: null,
        ],
    ]
    if (!method_options) {
        return analysis_rows.collectEntries { row -> [(row[0].id): defaults] }
    }

    def document_path = method_options.toString()
    def fail = { analysis_id, option, reason ->
        error("[nf-core/gwas] ERROR: Method-options document '${document_path}', analysis_id '${analysis_id}', option '${option}': ${reason}")
    }
    def document_file = file(method_options)
    if (!document_file.exists()) {
        fail.call('<document>', '<root>', 'file does not exist')
    }

    def document = null
    try {
        document = new groovy.json.JsonSlurper().parseText(document_file.text)
    }
    catch (Exception exception) {
        fail.call('<document>', '<root>', "malformed JSON (${exception.message})")
    }
    if (!(document instanceof Map)) {
        fail.call('<document>', '<root>', 'expected an object keyed by analysis_id')
    }

    def analyses = analysis_rows.collectEntries { row ->
        def meta = row[0]
        [(meta.id): [
            association_methods: tokenizeMethodSelector(meta.association_methods),
            heritability_methods: tokenizeMethodSelector(meta.heritability_methods),
            is_binary: meta.trait_type == 'binary',
        ]]
    }
    def resolved = analyses.collectEntries { analysis_id, _methods -> [(analysis_id): defaults] }

    document.each { analysis_id, families ->
        if (!analyses.containsKey(analysis_id)) {
            fail.call(analysis_id, '<analysis>', 'analysis identifier is not declared in the analysis manifest')
        }
        if (!(families instanceof Map)) {
            fail.call(analysis_id, '<analysis>', 'expected a method-family object')
        }
        def unknown_families = families.keySet().findAll { family -> !(family in ['gcta', 'ldak', 'regenie']) }
        if (unknown_families) {
            fail.call(analysis_id, unknown_families.first().toString(), 'unknown method family; accepted families are gcta, ldak and regenie')
        }

        def gcta = families.containsKey('gcta') ? families.gcta : [:]
        def ldak = families.containsKey('ldak') ? families.ldak : [:]
        def regenie = families.containsKey('regenie') ? families.regenie : [:]
        if (!(gcta instanceof Map)) {
            fail.call(analysis_id, 'gcta', 'expected an option object')
        }
        if (!(ldak instanceof Map)) {
            fail.call(analysis_id, 'ldak', 'expected an option object')
        }
        if (!(regenie instanceof Map)) {
            fail.call(analysis_id, 'regenie', 'expected an option object')
        }
        def accepted_regenie = ['step1_bsize', 'firth', 'firth_approx', 'firth_p_threshold', 'min_mac']
        def unknown_regenie = regenie.keySet().findAll { option -> !(option in accepted_regenie) }
        if (unknown_regenie) {
            def option = unknown_regenie.first()
            def reason = option in ['step2_bsize', 'step1_mode', 'step1_jobs', 'lowmem']
                ? 'operational tuning must be supplied through run/profile configuration'
                : "unknown option; accepted REGENIE options are ${accepted_regenie.join(', ')}"
            fail.call(analysis_id, "regenie.${option}", reason)
        }
        def accepted_ldak = [
            'model',
            'power',
            'weights_policy',
            'weights',
            'relatedness_filter',
            'kvik_step1_subset',
            'predictor_extract',
        ]
        def unknown_ldak = ldak.keySet().findAll { option -> !(option in accepted_ldak) }
        if (unknown_ldak) {
            def option = unknown_ldak.first()
            def reason = option in ['threads', 'jobs', 'partitions']
                ? 'operational tuning must be supplied through run/profile configuration'
                : "unknown option; accepted LDAK options are ${accepted_ldak.join(', ')}"
            fail.call(analysis_id, "ldak.${option}", reason)
        }

        def accepted_gcta = [
            'grm_maf',
            'grm_extract',
            'reml_no_constrain',
            'sparse_cutoff',
            'ld_score_region_kb',
            'ld_bins',
            'ldms_maf_edges',
        ]
        def unknown_gcta = gcta.keySet().findAll { option -> !(option in accepted_gcta) }
        if (unknown_gcta) {
            def option = unknown_gcta.first()
            def reason = option == 'gcta_grm_parts'
                ? 'partition count is operational and must be supplied through run/profile configuration'
                : "unknown option; accepted GCTA options are ${accepted_gcta.join(', ')}"
            fail.call(analysis_id, "gcta.${option}", reason)
        }

        def methods = analyses[analysis_id]
        def selects_regenie = methods.association_methods.contains('regenie')
        if (regenie && !selects_regenie) {
            fail.call(analysis_id, "regenie.${regenie.keySet().first()}", "analysis does not select 'regenie'")
        }
        def step1_bsize = regenie.containsKey('step1_bsize') ? regenie.step1_bsize : 1000
        if (!(step1_bsize instanceof Number) || step1_bsize < 1 || step1_bsize != step1_bsize.toInteger()) {
            fail.call(analysis_id, 'regenie.step1_bsize', 'expected one positive integer')
        }
        def firth = regenie.containsKey('firth') ? regenie.firth : true
        def firth_approx = regenie.containsKey('firth_approx') ? regenie.firth_approx : true
        def firth_p_threshold = regenie.containsKey('firth_p_threshold') ? regenie.firth_p_threshold : 0.01
        def min_mac = regenie.containsKey('min_mac') ? regenie.min_mac : null
        if (!(firth instanceof Boolean)) {
            fail.call(analysis_id, 'regenie.firth', 'expected a boolean')
        }
        if (!(firth_approx instanceof Boolean)) {
            fail.call(analysis_id, 'regenie.firth_approx', 'expected a boolean')
        }
        if (regenie.containsKey('firth_approx') && firth_approx && !firth) {
            fail.call(analysis_id, 'regenie.firth_approx', "requires 'regenie.firth' to be true")
        }
        if (!(firth_p_threshold instanceof Number) || firth_p_threshold <= 0 || firth_p_threshold > 1) {
            fail.call(analysis_id, 'regenie.firth_p_threshold', 'expected a number greater than 0 and at most 1')
        }
        if (min_mac != null && (!(min_mac instanceof Number) || min_mac < 0)) {
            fail.call(analysis_id, 'regenie.min_mac', 'expected a non-negative number or null')
        }
        ['firth', 'firth_approx', 'firth_p_threshold'].each { option ->
            if (regenie.containsKey(option) && !methods.is_binary) {
                fail.call(analysis_id, "regenie.${option}", 'option is consumed by binary-trait REGENIE analyses only')
            }
        }
        if (regenie.containsKey('firth_p_threshold') && !firth) {
            fail.call(analysis_id, 'regenie.firth_p_threshold', "requires 'regenie.firth' to be true")
        }
        def selects_gcta = methods.association_methods.contains('gcta_fastgwa') || methods.heritability_methods.any { method -> method in ['gcta_greml', 'gcta_greml_ldms'] }
        if (gcta && !selects_gcta) {
            fail.call(analysis_id, "gcta.${gcta.keySet().first()}", 'analysis does not select a GCTA method')
        }
        if (gcta.containsKey('reml_no_constrain') && !methods.heritability_methods.any { method -> method in ['gcta_greml', 'gcta_greml_ldms'] }) {
            fail.call(analysis_id, 'gcta.reml_no_constrain', "option is consumed by GCTA GREML estimators only, but this analysis selects neither 'gcta_greml' nor 'gcta_greml_ldms'")
        }
        ['grm_maf', 'grm_extract'].each { option ->
            if (gcta.containsKey(option) && !methods.heritability_methods.contains('gcta_greml')) {
                fail.call(analysis_id, "gcta.${option}", "option is consumed by 'gcta_greml' only, which this analysis does not select")
            }
        }
        if (gcta.containsKey('sparse_cutoff') && !methods.association_methods.contains('gcta_fastgwa')) {
            fail.call(analysis_id, 'gcta.sparse_cutoff', "option is consumed by 'gcta_fastgwa' only, which this analysis does not select")
        }
        ['ld_score_region_kb', 'ld_bins', 'ldms_maf_edges'].each { option ->
            if (gcta.containsKey(option) && !methods.heritability_methods.contains('gcta_greml_ldms')) {
                fail.call(analysis_id, "gcta.${option}", "option is consumed by 'gcta_greml_ldms' only, which this analysis does not select")
            }
        }

        def grm_maf = gcta.containsKey('grm_maf') ? gcta.grm_maf : null
        if (grm_maf != null && (!(grm_maf instanceof Number) || grm_maf < 0 || grm_maf > 0.5)) {
            fail.call(analysis_id, 'gcta.grm_maf', 'expected a number between 0 and 0.5 inclusive')
        }
        def sparse_cutoff = gcta.containsKey('sparse_cutoff') ? gcta.sparse_cutoff : 0.05
        if (!(sparse_cutoff instanceof Number) || sparse_cutoff < 0 || sparse_cutoff > 1) {
            fail.call(analysis_id, 'gcta.sparse_cutoff', 'expected a number between 0 and 1 inclusive')
        }
        def ld_score_region_kb = gcta.containsKey('ld_score_region_kb') ? gcta.ld_score_region_kb : 200
        if (!(ld_score_region_kb instanceof Number) || ld_score_region_kb < 1 || ld_score_region_kb != ld_score_region_kb.toInteger()) {
            fail.call(analysis_id, 'gcta.ld_score_region_kb', 'expected one positive integer')
        }
        def ld_bins = gcta.containsKey('ld_bins') ? gcta.ld_bins : 4
        if (!(ld_bins instanceof Number) || ld_bins < 1 || ld_bins != ld_bins.toInteger()) {
            fail.call(analysis_id, 'gcta.ld_bins', 'expected one positive integer')
        }
        def ldms_maf_edges = gcta.containsKey('ldms_maf_edges') ? gcta.ldms_maf_edges : [0, 0.01, 0.05, 0.2, 0.5]
        if (!(ldms_maf_edges instanceof List) || ldms_maf_edges.size() < 2 || !ldms_maf_edges.every { edge -> edge instanceof Number && edge >= 0 && edge <= 0.5 }) {
            fail.call(analysis_id, 'gcta.ldms_maf_edges', 'expected a numeric boundary list spanning 0 to 0.5')
        }
        if (ldms_maf_edges.first() != 0 || ldms_maf_edges.last() != 0.5) {
            fail.call(analysis_id, 'gcta.ldms_maf_edges', 'boundaries must start at 0 and end at 0.5')
        }
        if ((1..<ldms_maf_edges.size()).any { index -> ldms_maf_edges[index] <= ldms_maf_edges[index - 1] }) {
            fail.call(analysis_id, 'gcta.ldms_maf_edges', 'boundaries must be strictly increasing')
        }
        def reml_no_constrain = gcta.containsKey('reml_no_constrain') ? gcta.reml_no_constrain : false
        if (!(reml_no_constrain instanceof Boolean)) {
            fail.call(analysis_id, 'gcta.reml_no_constrain', 'expected a boolean')
        }
        def grm_extract = []
        if (gcta.containsKey('grm_extract')) {
            if (!(gcta.grm_extract instanceof String) || !gcta.grm_extract.trim()) {
                fail.call(analysis_id, 'gcta.grm_extract', 'expected a non-empty resource path string')
            }
            grm_extract = file(gcta.grm_extract)
            if (!grm_extract.exists()) {
                fail.call(analysis_id, 'gcta.grm_extract', "resource path '${gcta.grm_extract}' does not exist")
            }
        }

        def model = ldak.containsKey('model') ? ldak.model : 'human_default'
        if (!(model in ['human_default', 'custom'])) {
            fail.call(analysis_id, 'ldak.model', "expected 'human_default' or 'custom'")
        }
        def power = ldak.containsKey('power') ? ldak.power : -0.25
        if (!(power instanceof Number) || power < -2 || power > 0) {
            fail.call(analysis_id, 'ldak.power', 'expected a number between -2 and 0 inclusive')
        }
        if (model == 'human_default' && power != -0.25) {
            fail.call(analysis_id, 'ldak.power', "model 'human_default' fixes power at -0.25; set model to 'custom' to supply another power")
        }
        def weights_policy = ldak.containsKey('weights_policy') ? ldak.weights_policy : 'equal'
        if (!(weights_policy in ['equal', 'default', 'provided'])) {
            fail.call(analysis_id, 'ldak.weights_policy', "expected 'equal', 'default' or 'provided'")
        }
        if (ldak.containsKey('weights') && weights_policy != 'provided') {
            fail.call(analysis_id, 'ldak.weights', "resource is only accepted when weights_policy is 'provided', got '${weights_policy}'")
        }
        if (weights_policy == 'provided' && !ldak.containsKey('weights')) {
            fail.call(analysis_id, 'ldak.weights', "weights_policy is 'provided' but no resource is supplied")
        }
        def weights = []
        if (ldak.containsKey('weights')) {
            if (!(ldak.weights instanceof String) || !ldak.weights.trim()) {
                fail.call(analysis_id, 'ldak.weights', 'expected a non-empty resource path string')
            }
            weights = file(ldak.weights)
            if (!weights.exists()) {
                fail.call(analysis_id, 'ldak.weights', "resource path '${ldak.weights}' does not exist")
            }
        }
        def relatedness_filter = ldak.containsKey('relatedness_filter') ? ldak.relatedness_filter : false
        if (!(relatedness_filter instanceof Boolean)) {
            fail.call(analysis_id, 'ldak.relatedness_filter', 'expected a boolean')
        }
        def kvik_step1_subset = ldak.containsKey('kvik_step1_subset') ? ldak.kvik_step1_subset : 'all'
        if (!(kvik_step1_subset in ['all', 'thin_common', 'provided'])) {
            fail.call(analysis_id, 'ldak.kvik_step1_subset', "expected 'all', 'thin_common' or 'provided'")
        }
        if (ldak.containsKey('predictor_extract') && kvik_step1_subset != 'provided') {
            fail.call(analysis_id, 'ldak.predictor_extract', "resource is only accepted when kvik_step1_subset is 'provided', got '${kvik_step1_subset}'")
        }
        if (kvik_step1_subset == 'provided' && !ldak.containsKey('predictor_extract')) {
            fail.call(analysis_id, 'ldak.predictor_extract', "kvik_step1_subset is 'provided' but no resource is supplied")
        }
        def predictor_extract = []
        if (ldak.containsKey('predictor_extract')) {
            if (!(ldak.predictor_extract instanceof String) || !ldak.predictor_extract.trim()) {
                fail.call(analysis_id, 'ldak.predictor_extract', 'expected a non-empty resource path string')
            }
            predictor_extract = file(ldak.predictor_extract)
            if (!predictor_extract.exists()) {
                fail.call(analysis_id, 'ldak.predictor_extract', "resource path '${ldak.predictor_extract}' does not exist")
            }
        }
        def selects_ldak_kinship = methods.heritability_methods.any { method -> method in ['ldak_reml', 'ldak_he', 'ldak_pcgc'] }
        def selects_ldak_kvik = methods.association_methods.contains('ldak_kvik')
        if (ldak && !selects_ldak_kinship && !selects_ldak_kvik) {
            fail.call(analysis_id, "ldak.${ldak.keySet().first()}", 'analysis does not select an LDAK method')
        }
        if (ldak.containsKey('weights') && !selects_ldak_kinship) {
            fail.call(analysis_id, 'ldak.weights', 'resource is consumed by LDAK kinship methods only, which this analysis does not select')
        }
        ['model', 'power', 'weights_policy', 'relatedness_filter'].each { option ->
            if (ldak.containsKey(option) && !selects_ldak_kinship) {
                fail.call(analysis_id, "ldak.${option}", 'option is consumed by LDAK kinship methods only, which this analysis does not select')
            }
        }
        ['kvik_step1_subset', 'predictor_extract'].each { option ->
            if (ldak.containsKey(option) && !selects_ldak_kvik) {
                fail.call(analysis_id, "ldak.${option}", "option is consumed by 'ldak_kvik' only, which this analysis does not select")
            }
        }


        resolved[analysis_id] = [
            gcta: [
                grm_maf: grm_maf,
                grm_extract: grm_extract,
                reml_no_constrain: reml_no_constrain,
                sparse_cutoff: sparse_cutoff,
                ld_score_region_kb: ld_score_region_kb,
                ld_bins: ld_bins,
                ldms_maf_edges: ldms_maf_edges,
            ],
            ldak: [
                model: model,
                power: power,
                weights_policy: weights_policy,
                weights: weights,
                relatedness_filter: relatedness_filter,
                kvik_step1_subset: kvik_step1_subset,
                predictor_extract: predictor_extract,
            ],
            regenie: [
                step1_bsize: step1_bsize.toInteger(),
                firth: firth,
                firth_approx: firth_approx,
                firth_p_threshold: firth_p_threshold,
                min_mac: min_mac,
            ],
        ]
    }
    return resolved
}

//
// Validate the linked manifests as whole lists, join cohort-owned facts to every analysis through
// cohort_id, and build the canonical tuple consumed by GWAS. The maps retain the first occurrence
// only long enough to diagnose later duplicates; any error aborts before the list becomes a channel.
//
def validateRelationalInput(cohort_rows, analysis_rows, cohort_manifest, analysis_manifest, cohort_schema, analysis_schema, method_options = null) {
    def cohort_columns = getSamplesheetPositionalColumns(cohort_schema)
    def analysis_columns = getSamplesheetPositionalColumns(analysis_schema)

    def errors = []
    def cohorts_by_id = [:]
    def line_by_analysis_id = [:]
    def validated_rows = []
    def method_options_by_analysis = validateMethodOptions(method_options, analysis_rows)

    cohort_rows.eachWithIndex { row, index ->
        def line = index + 2
        def cohort_meta = row[0]
        def cells = [cohort_columns, row[1..-1]].transpose().collectEntries()
        def cohort_id = cohort_meta.cohort
        def reject = { field, message ->
            def named = field instanceof List ? field : [field]
            def label = named.size() > 1
                ? "fields ${named.collect { name -> "'${name}'" }.join(', ')}"
                : "field '${named.first()}'"
            errors << "  - ${cohort_manifest} row ${line} (cohort_id '${cohort_id}'), ${label}: ${message}"
        }

        def genotype_format = validateGenotypeGroup(cells, reject)
        def genotype_files = genotype_format
            ? getGenotypeGroups()[genotype_format].collect { column -> cells[column] }
            : []
        def definition = [
            genome_build: cohort_meta.build,
            ancestry: cohort_meta.ancestry,
        ] + getGenotypeGroups().values().flatten().collectEntries { field -> [(field): normaliseCellValue(cells[field])?.toString() ?: ''] }
        def known = cohorts_by_id[cohort_id]

        if (known) {
            if (known.definition == definition) {
                reject.call('cohort_id', "duplicate cohort_id '${cohort_id}', identical definition on row ${known.line}")
            }
            else {
                def differing_fields = definition.keySet().findAll { field -> known.definition[field] != definition[field] }
                reject.call('cohort_id', "duplicate cohort_id '${cohort_id}' conflicts with row ${known.line} in field${differing_fields.size() > 1 ? 's' : ''} ${differing_fields.collect { field -> "'${field}'" }.join(', ')}")
            }
        }
        else {
            cohorts_by_id[cohort_id] = [
                line: line,
                meta: cohort_meta,
                genotype_format: genotype_format,
                genotype_files: genotype_files,
                definition: definition,
            ]
        }
    }

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
        validateMethodSelectors(association_methods, heritability_methods, reject)

        def routes = getMethodRoutes(association_methods, heritability_methods)
        def settings = getAnalysisSettings(analysis_meta)
        def is_binary = analysis_meta.trait_type == 'binary'
        validateTraitColumns(is_binary, settings, reject)
        validateMethodConditionedColumns(settings, routes, reject)

        if (cohort) {
            validated_rows << [
                analysis_meta + [
                    build: cohort.meta.build,
                    ancestry: cohort.meta.ancestry,
                    association_methods: association_methods,
                    heritability_methods: heritability_methods,
                    genotype_format: cohort.genotype_format,
                    is_binary: is_binary,
                    has_covariates: cells.quant_covariates || cells.cat_covariates ? true : false,
                    case_value: settings.case_value == null ? null : settings.case_value.toString(),
                    control_value: settings.control_value == null ? null : settings.control_value.toString(),
                    population_prevalence: settings.population_prevalence,
                    method_options: method_options_by_analysis[analysis_id],
                ],
                cohort.genotype_files,
                cells.phenotype,
                cells.quant_covariates,
                cells.cat_covariates,
                method_options_by_analysis[analysis_id].ldak.predictor_extract,
                method_options_by_analysis[analysis_id].ldak.weights,
            ]
        }
    }

    if (errors) {
        error("[nf-core/gwas] ERROR: Validation of linked manifests failed!\n\n${errors.join('\n')}\n")
    }

    return validated_rows
}


//
// Build the MultiQC analysis-plan table from validated, joined manifest metadata.
//

def analysisPlanJson(mqc_analysis_plan_yaml, analysis_metadata) {
    def static_config = new org.yaml.snakeyaml.Yaml().load(mqc_analysis_plan_yaml.text)
    def header_keys = static_config.headers.keySet()
    def data = analysis_metadata
        .toSorted { meta -> meta.id }
        .collectEntries { meta ->
            def raw = [
                cohort: meta.cohort,
                trait: meta.trait,
                trait_type: meta.trait_type,
                genome_build: meta.build,
                ancestry: meta.ancestry,
                association_methods: meta.association_methods ? meta.association_methods.join(', ') : '-',
                heritability_methods: meta.heritability_methods ? meta.heritability_methods.join(', ') : '-',
            ]
            def undeclared = raw.keySet() - header_keys
            def unpopulated = header_keys - raw.keySet()
            if (undeclared || unpopulated) {
                error("multiqc_analysis_plan.yml contract mismatch; generated columns missing from headers: ${undeclared}; declared headers without generated data: ${unpopulated}")
            }

            def cells = [:]
            header_keys.each { key -> cells[key] = raw[key] }
            [(meta.id): cells]
        }

    return groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(static_config + [data: data]))
}


//
// Generate a route-aware methods description for MultiQC.
//

def methodsDescriptionText(mqc_methods_yaml, selected_methods = [association: [], heritability: []]) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>No version-specific pipeline DOI was declared for this build.</li>"

    // Tool references
    meta["tool_citations"] = toolCitationText(selected_methods)
    meta["tool_bibliography"] = toolBibliographyText(selected_methods)
    meta["command_line"] = escapeHtml(workflow.commandLine)

    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}


def selectedCitationKeys(selected_methods) {
    def association = (selected_methods.association ?: []) as Set
    def heritability = (selected_methods.heritability ?: []) as Set
    def known_association = ['plink2', 'regenie', 'gcta_fastgwa', 'ldak_kvik'] as Set
    def known_heritability = ['gcta_greml', 'gcta_greml_ldms', 'ldak_reml', 'ldak_he', 'ldak_pcgc'] as Set
    def unknown = (association - known_association) + (heritability - known_heritability)
    if (unknown) {
        error("Cannot generate methods citations for unknown method selectors: ${unknown.toList().sort().join(', ')}")
    }

    def keys = []
    if ('plink2' in association) {
        keys << 'plink2'
    }
    if ('regenie' in association) {
        keys << 'regenie'
    }
    if ('gcta_fastgwa' in association) {
        keys << 'gcta_fastgwa'
    }
    if ('gcta_greml' in heritability) {
        keys << 'gcta_greml'
    }
    if ('gcta_greml_ldms' in heritability) {
        keys << 'gcta_greml_ldms'
    }
    if ('ldak_kvik' in association) {
        keys << 'ldak_kvik'
    }
    if (heritability.intersect(['ldak_reml', 'ldak_he', 'ldak_pcgc'])) {
        keys << 'ldak'
    }
    if (association) {
        keys << 'gwaslab'
    }
    keys << 'multiqc'
    return keys
}


def toolCitationText(selected_methods) {
    def keys = selectedCitationKeys(selected_methods)
    def association_labels = [
        plink2: 'PLINK 2 (Chang <em>et al.</em>, 2015)',
        regenie: 'REGENIE (Mbatchou <em>et al.</em>, 2021)',
        gcta_fastgwa: 'GCTA fastGWA (Jiang <em>et al.</em>, 2019)',
        ldak_kvik: 'LDAK-KVIK (Hof and Speed, 2025)',
    ]
    def heritability_labels = [
        gcta_greml: 'GCTA GREML (Yang <em>et al.</em>, 2011)',
        gcta_greml_ldms: 'GCTA GREML-LDMS (Yang <em>et al.</em>, 2015)',
        ldak: 'LDAK (Speed <em>et al.</em>, 2012)',
    ]
    def sentences = []
    def selected_association = keys.findAll { key -> association_labels.containsKey(key) }.collect { key -> association_labels[key] }
    def selected_heritability = keys.findAll { key -> heritability_labels.containsKey(key) }.collect { key -> heritability_labels[key] }
    if (selected_association) {
        sentences << "Association testing was performed with ${joinProseList(selected_association)}."
        sentences << "Association summary statistics were harmonised with GWASLab."
    }
    if (selected_heritability) {
        sentences << "SNP-based heritability was estimated with ${joinProseList(selected_heritability)}."
    }
    sentences << "The run report was generated with MultiQC (Ewels <em>et al.</em>, 2016)."
    return sentences.join(' ')
}


def toolBibliographyText(selected_methods) {
    def bibliography = [
        plink2: '<li>Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. <em>GigaScience</em>. 2015;4:7. doi: <a href="https://doi.org/10.1186/s13742-015-0047-8">10.1186/s13742-015-0047-8</a>.</li>',
        regenie: '<li>Mbatchou J, Barnard L, Backman J, et al. Computationally efficient whole-genome regression for quantitative and binary traits. <em>Nature Genetics</em>. 2021;53:1097-1103. doi: <a href="https://doi.org/10.1038/s41588-021-00870-7">10.1038/s41588-021-00870-7</a>.</li>',
        gcta_fastgwa: '<li>Jiang L, Zheng Z, Qi T, et al. A resource-efficient tool for mixed model association analysis of large-scale data. <em>Nature Genetics</em>. 2019;51:1749-1755. doi: <a href="https://doi.org/10.1038/s41588-019-0530-8">10.1038/s41588-019-0530-8</a>.</li>',
        gcta_greml: '<li>Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: a tool for genome-wide complex trait analysis. <em>American Journal of Human Genetics</em>. 2011;88:76-82. doi: <a href="https://doi.org/10.1016/j.ajhg.2010.11.011">10.1016/j.ajhg.2010.11.011</a>.</li>',
        gcta_greml_ldms: '<li>Yang J, Bakshi A, Zhu Z, et al. Genetic variance estimation with imputed variants finds negligible missing heritability for human height and body mass index. <em>Nature Genetics</em>. 2015;47:1114-1120. doi: <a href="https://doi.org/10.1038/ng.3390">10.1038/ng.3390</a>.</li>',
        ldak_kvik: '<li>Hof JP, Speed D. LDAK-KVIK performs fast and powerful mixed-model association analysis of quantitative and binary phenotypes. <em>Nature Genetics</em>. 2025;57:2116-2123. doi: <a href="https://doi.org/10.1038/s41588-025-02286-z">10.1038/s41588-025-02286-z</a>.</li>',
        ldak: '<li>Speed D, Hemani G, Johnson MR, Balding DJ. Improved heritability estimation from genome-wide SNPs. <em>American Journal of Human Genetics</em>. 2012;91:1011-1021. doi: <a href="https://doi.org/10.1016/j.ajhg.2012.10.010">10.1016/j.ajhg.2012.10.010</a>.</li>',
        gwaslab: '<li>GWASLab. <a href="https://cloufield.github.io/gwaslab/">https://cloufield.github.io/gwaslab/</a>.</li>',
        multiqc: '<li>Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. <em>Bioinformatics</em>. 2016;32:3047-3048. doi: <a href="https://doi.org/10.1093/bioinformatics/btw354">10.1093/bioinformatics/btw354</a>.</li>',
    ]
    return selectedCitationKeys(selected_methods).collect { key -> bibliography[key] }.join('\n    ')
}


def joinProseList(items) {
    if (items.size() == 1) {
        return items[0]
    }
    if (items.size() == 2) {
        return items.join(' and ')
    }
    return "${items[0..-2].join(', ')}, and ${items[-1]}"
}


def escapeHtml(value) {
    return value
        .toString()
        .replace('&', '&amp;')
        .replace('<', '&lt;')
        .replace('>', '&gt;')
        .replace('"', '&quot;')
        .replace("'", '&#39;')
}


// Any file-like input contributing to scientific identity is hashed by content. Byte-identical files
// materialised for different analysis IDs are valid reuse candidates; equal basenames are not sufficient.
def digestScientificInput(input_file) {
    if (!input_file) {
        return 'absent'
    }
    def input_path = input_file instanceof java.nio.file.Path ? input_file : input_file.toPath()
    def digest = java.security.MessageDigest.getInstance('SHA-256')
    java.nio.file.Files
        .newInputStream(input_path)
        .withCloseable { input ->
            input.eachByte(8192) { buffer, count ->
                digest.update(buffer, 0, count)
            }
        }
    return digest.digest().encodeHex().toString()
}


// Prediction reuse keys share one deterministic map serialisation and the same 12-character SHA-256 prefix.
def buildCanonicalPredictionKey(identity) {
    def canonical = identity
        .sort { entry -> entry.key }
        .collect { name, value -> "${name}=${value}" }
        .join('\n')
    return java.security.MessageDigest
        .getInstance('SHA-256')
        .digest(canonical.getBytes('UTF-8'))
        .encodeHex()
        .toString()
        .substring(0, 12)
}


// REGENIE Step 1 reuse requires every cohort-defining and scientific input to agree. Execution-only
// controls are intentionally absent; Step 1 block size remains because it changes the fitted model.
def buildRegeniePredictionKey(meta, phenotype, covariates, step1_bsize) {
    def identity = [
        cohort: meta.cohort,
        trait: meta.trait,
        is_binary: meta.is_binary,
        phenotype: digestScientificInput(phenotype),
        covariates: digestScientificInput(covariates),
        step1_bsize: step1_bsize,
    ]
    return buildCanonicalPredictionKey(identity)
}


// LDAK-KVIK Step 1 reuse additionally depends on predictor policy and optional predictor-list bytes.
def buildKvikPredictionKey(meta, phenotype, quant_covariates, cat_covariates, subset_policy, predictor_extract) {
    def identity = [
        cohort: meta.cohort,
        trait: meta.trait,
        is_binary: meta.is_binary,
        phenotype: digestScientificInput(phenotype),
        quant_covariates: digestScientificInput(quant_covariates),
        cat_covariates: digestScientificInput(cat_covariates),
        subset_policy: subset_policy,
        predictor_extract: digestScientificInput(predictor_extract),
    ]
    return buildCanonicalPredictionKey(identity)
}
