// Pipeline-specific initialisation, completion, reporting and shared identity primitives.
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

// SUBWORKFLOW: Local to the pipeline
include { VALIDATE_GWAS_INPUT     } from '../validate_gwas_input'

// FUNCTION: Local to the pipeline
include { getMethodCapabilities   } from '../validate_gwas_input/method_registry'

// PLUGIN
include { paramsSummaryMap        } from 'plugin/nf-schema'

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
    // SUBWORKFLOW: Validate linked input manifests and construct canonical analysis records
    //
    VALIDATE_GWAS_INPUT(
        cohort_manifest,
        analysis_manifest,
        method_options,
    )

    emit:
    analyses = VALIDATE_GWAS_INPUT.out.analyses // channel: [ val(meta), [ path(genotype_file), ... ], path(phenotype), path(quant_covariates), path(cat_covariates) ]
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
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta['manifest_map'] = workflow.manifest.toMap()

    if (meta.manifest_map.doi) {
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(',')
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace('https://doi.org/', '').replace(' ', '')}\'>${doi_ref.replace('https://doi.org/', '').replace(' ', '')}</a>), "
        }
        meta['doi_text'] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta['doi_text'] = ""
    }
    meta['nodoi_text'] = meta.manifest_map.doi ? "" : '<li>No version-specific pipeline DOI was declared for this build.</li>'
    meta['tool_citations'] = toolCitationText(selected_methods)
    meta['tool_bibliography'] = toolBibliographyText(selected_methods)
    meta['command_line'] = escapeHtml(workflow.commandLine)

    def methods_text = mqc_methods_yaml.text
    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)
    return description_html.toString()
}

// The citation keys a run's selected methods carry, in the registry's declared order, followed by the
// report renderer. Every selected token must be registered: a token the registry does not know has no
// citation to render, and an unknown token at this point is a routing error rather than a user error.
def selectedCitationKeys(selected_methods) {
    def association = (selected_methods.association ?: []) as Set
    def heritability = (selected_methods.heritability ?: []) as Set
    def capabilities = getMethodCapabilities()
    def known_association = capabilities.findAll { _token, details -> details.domain == 'association' }.keySet() as Set
    def known_heritability = capabilities.findAll { _token, details -> details.domain == 'heritability' }.keySet() as Set
    def unknown = (association - known_association) + (heritability - known_heritability)
    if (unknown) {
        error("Cannot generate methods citations for unknown method selectors: ${unknown.toList().sort().join(', ')}")
    }

    def keys = (association + heritability)
        .collectMany { token -> capabilities[token].citation_keys }
        .findAll { key -> key }
        .unique()
    keys << 'multiqc'
    return keys
}

// Each wired route adds its prose label here, keyed by selector token, together with its bibliography entry
// below, keyed by the citation key its registry entry declares.
def getMethodCitationLabels() {
    return [
        association: [:],
        heritability: [:],
    ]
}

def toolCitationText(selected_methods) {
    def association = (selected_methods.association ?: []) as Set
    def heritability = (selected_methods.heritability ?: []) as Set
    def labels = getMethodCitationLabels()
    def sentences = []
    def selected_association = labels.association.findAll { token, _label -> token in association }.values().toList()
    def selected_heritability = labels.heritability.findAll { token, _label -> token in heritability }.values().toList().unique()
    if (selected_association) {
        sentences << "Association testing was performed with ${joinProseList(selected_association)}."
    }
    if (selected_heritability) {
        sentences << "SNP-based heritability was estimated with ${joinProseList(selected_heritability)}."
    }
    sentences << 'The run report was generated with MultiQC (Ewels <em>et al.</em>, 2016).'
    return sentences.join(' ')
}

def toolBibliographyText(selected_methods) {
    def bibliography = [
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

// Any present file-like input contributing to scientific identity is hashed by content. Byte-identical
// files materialised for different analysis IDs are valid reuse candidates; equal basenames are not sufficient.
def digestFileBytes(input_file) {
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

// Reuse keys share one short SHA-256 primitive over already-canonical UTF-8 text.
def digestIdentityText(canonical) {
    return java.security.MessageDigest
        .getInstance('SHA-256')
        .digest(canonical.getBytes('UTF-8'))
        .encodeHex()
        .toString()
        .substring(0, 12)
}

// Scientific values canonicalise numeric spellings while preserving identifiers exactly as declared.
def canonicaliseScientificValue(value) {
    if (value == null) {
        return 'null'
    }
    if (value instanceof Boolean) {
        return value ? 'true' : 'false'
    }
    if (value instanceof Map) {
        return '{' + value.sort { entry -> entry.key }.collect { name, entry -> "${name}=${canonicaliseScientificValue(entry)}" }.join(',') + '}'
    }
    if (value instanceof Collection) {
        return '[' + value.collect { entry -> canonicaliseScientificValue(entry) }.join(',') + ']'
    }
    if (value instanceof Number) {
        return new BigDecimal(value.toString()).stripTrailingZeros().toPlainString()
    }
    def text = value.toString().trim()
    return text ==~ /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/
        ? new BigDecimal(text).stripTrailingZeros().toPlainString()
        : text
}

def canonicaliseScientificIdentifier(value) {
    if (value == null) {
        return 'null'
    }
    if (value instanceof Map) {
        return '{' + value.sort { entry -> entry.key }.collect { name, entry -> "${name}=${canonicaliseScientificIdentifier(entry)}" }.join(',') + '}'
    }
    if (value instanceof Collection) {
        return '[' + value.collect { entry -> canonicaliseScientificIdentifier(entry) }.join(',') + ']'
    }
    return value.toString().trim()
}

// Artifact keys contain only immutable parent/input identities and effective scientific settings.
def buildScientificArtifactKey(identity, settings) {
    def rendered = identity.collectEntries { name, value -> [(name): canonicaliseScientificIdentifier(value)] } + [settings: canonicaliseScientificValue(settings)]
    def canonical = rendered
        .sort { entry -> entry.key }
        .collect { name, text -> "${name}=${text}" }
        .join('\n')
    return digestIdentityText(canonical)
}

// Prediction reuse keys share one deterministic map serialisation.
def buildCanonicalPredictionKey(identity) {
    def canonical = identity
        .sort { entry -> entry.key }
        .collect { name, value -> "${name}=${value}" }
        .join('\n')
    return digestIdentityText(canonical)
}
