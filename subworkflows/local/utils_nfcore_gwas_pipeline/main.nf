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
    summary_statistics_manifest // channel: val(summary_statistics_manifest)
    relationship_manifest // channel: val(relationship_manifest)
    reference_catalog // channel: val(reference_catalog)
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

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> [--cohort_manifest cohorts.csv --analysis_manifest analyses.csv | --summary_statistics_manifest summaries.csv] --outdir <OUTDIR>"

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
        summary_statistics_manifest,
        relationship_manifest,
        reference_catalog,
        method_options,
    )

    emit:
    analyses           = VALIDATE_GWAS_INPUT.out.analyses // channel: [ val(meta), [ path(genotype_file), ... ], path(phenotype), path(quant_covariates), path(cat_covariates), path(kvik_extract), path(ldak_weights) ]
    summary_statistics = VALIDATE_GWAS_INPUT.out.summary_statistics // channel: [ val(meta), path(source) ]
    relationships      = VALIDATE_GWAS_INPUT.out.relationships // channel: [ val(meta), [ path(genotype_file), ... ], path(pair_quant_covariates), path(pair_cat_covariates) ]
    unary_requests     = VALIDATE_GWAS_INPUT.out.unary_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ]
    meta_requests      = VALIDATE_GWAS_INPUT.out.meta_requests // channel: [ val(meta), val(source_summary_statistics_ids) ]
    pair_requests      = VALIDATE_GWAS_INPUT.out.pair_requests // channel: [ val(meta), path(hapmap3_snplist), path(reference_ld_scores), path(regression_weights), path(tagging_file) ]
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
// The optional GWASLab reference resources, keyed on genome build and on nothing else. The samplesheet's
// `genome_build` column is the only real heterogeneity in this set, so the reference selection reduces to
// this map. A build with no configured resource yields empty lists and harmonisation proceeds without one.
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
def methodsDescriptionText(mqc_methods_yaml, selected_methods = [association: [], heritability: [], pairwise: []]) {
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

def selectedCitationKeys(selected_methods) {
    def association = (selected_methods.association ?: []) as Set
    def heritability = (selected_methods.heritability ?: []) as Set
    def pairwise = (selected_methods.pairwise ?: []) as Set
    def summary_set = (selected_methods.summary_set ?: []) as Set
    def capabilities = getMethodCapabilities()
    def known_association = capabilities.findAll { _token, details -> details.domain == 'association' }.keySet() as Set
    def known_heritability = capabilities.findAll { _token, details -> details.domain in ['heritability', 'summary_unary'] }.keySet() as Set
    def known_pairwise = capabilities.findAll { _token, details -> details.domain == 'pairwise' }.keySet() as Set
    def known_summary_set = capabilities.findAll { _token, details -> details.domain == 'summary_set' }.keySet() as Set
    def unknown = (association - known_association) + (heritability - known_heritability) + (pairwise - known_pairwise) + (summary_set - known_summary_set)
    if (unknown) {
        error("Cannot generate methods citations for unknown method selectors: ${unknown.toList().sort().join(', ')}")
    }

    def citation_order = ['regenie', 'gcta_fastgwa', 'gcta_greml', 'gcta_greml_ldms', 'gcta_bivariate_reml', 'gcta_hereg', 'ldak_kvik', 'ldak', 'rhe_mc', 'mph', 'ldak_sumstats', 'ldsc', 'gwaslab', 'metasoft', 'mrmega']
    def keys = (association + heritability + pairwise + summary_set)
        .collectMany { token -> capabilities[token].citation_keys }
        .findAll { key -> key }
        .unique()
        .sort { key -> citation_order.indexOf(key) }
    if (association && !keys.contains('gwaslab')) {
        keys << 'gwaslab'
    }
    keys << 'multiqc'
    return keys
}

def toolCitationText(selected_methods) {
    def association = (selected_methods.association ?: []) as Set
    def heritability = (selected_methods.heritability ?: []) as Set
    def pairwise = (selected_methods.pairwise ?: []) as Set
    def summary_set = (selected_methods.summary_set ?: []) as Set
    def association_labels = [
        regenie: 'REGENIE (Mbatchou <em>et al.</em>, 2021)',
        gcta_fastgwa: 'GCTA fastGWA (Jiang <em>et al.</em>, 2019)',
        ldak_kvik: 'LDAK-KVIK (Hof and Speed, 2025)',
    ]
    def heritability_labels = [
        gcta_greml: 'GCTA GREML (Yang <em>et al.</em>, 2011)',
        gcta_greml_ldms: 'GCTA GREML-LDMS (Yang <em>et al.</em>, 2015)',
        ldak_reml: 'LDAK REML (Speed <em>et al.</em>, 2012)',
        ldak_he: 'LDAK Haseman-Elston regression (Speed <em>et al.</em>, 2012)',
        ldak_pcgc: 'LDAK PCGC regression (Speed <em>et al.</em>, 2012)',
        ldak_fast_he: 'LDAK fast Haseman-Elston regression (Speed <em>et al.</em>, 2012; Pazokitoroudi <em>et al.</em>, 2020)',
        ldak_fast_pcgc: 'LDAK fast PCGC regression (Speed <em>et al.</em>, 2012; Pazokitoroudi <em>et al.</em>, 2020)',
        mph_reml: 'MPH REML (Jiang, 2024)',
        mph_reml_ldms: 'MPH REML over LD- and MAF-stratified components (Jiang, 2024; Yang <em>et al.</em>, 2015)',
        ldak_sumher: 'LDAK SumHer (Speed and Balding, 2019)',
        ldsc_h2: 'LDSC (Bulik-Sullivan <em>et al.</em>, 2015)',
    ]
    def pairwise_labels = [
        gcta_bivariate_reml: 'GCTA bivariate REML (Lee <em>et al.</em>, 2012)',
        gcta_bivariate_reml_ldms: 'GCTA bivariate REML-LDMS (Lee <em>et al.</em>, 2012; Yang <em>et al.</em>, 2015)',
        gcta_bivariate_he: 'GCTA bivariate Haseman-Elston regression (Chen, 2014)',
        gcta_bivariate_he_ldms: 'GCTA bivariate Haseman-Elston regression over LD- and MAF-stratified components (Chen, 2014; Yang <em>et al.</em>, 2015)',
        mph_bivariate_reml: 'MPH bivariate REML (Jiang, 2024)',
        mph_bivariate_reml_ldms: 'MPH bivariate REML over LD- and MAF-stratified components (Jiang, 2024; Yang <em>et al.</em>, 2015)',
        ldak_sumcors: 'LDAK SumCors (Speed and Balding, 2019)',
        ldsc_rg: 'LDSC genetic correlation (Bulik-Sullivan <em>et al.</em>, 2015)',
    ]
    def meta_analysis_labels = [
        fixed: 'GWASLab inverse-variance fixed effects',
        random: 'GWASLab DerSimonian–Laird random effects',
        re2: 'METASOFT Han–Eskin RE2 (Han and Eskin, 2011)',
        mrmega: 'MR-MEGA (Mägi <em>et al.</em>, 2017)',
    ]
    def sentences = []
    def selected_association = association_labels.findAll { token, _label -> token in association }.values().toList()
    def selected_heritability = heritability_labels.findAll { token, _label -> token in heritability }.values().toList().unique()
    def selected_pairwise = pairwise_labels.findAll { token, _label -> token in pairwise }.values().toList()
    if (selected_association) {
        sentences << "Association testing was performed with ${joinProseList(selected_association)}."
        sentences << 'Association summary statistics were harmonised with GWASLab.'
    }
    if (selected_heritability) {
        sentences << "SNP-based heritability was estimated with ${joinProseList(selected_heritability)}."
    }
    if (selected_pairwise) {
        sentences << "Pairwise genetic covariance and correlation were estimated with ${joinProseList(selected_pairwise)}."
    }
    def selected_meta_analysis = meta_analysis_labels.findAll { token, _label -> token in summary_set }.values().toList()
    if (selected_meta_analysis) {
        sentences << "Common-variant meta-analysis used ${joinProseList(selected_meta_analysis)}."
    }
    sentences << 'The run report was generated with MultiQC (Ewels <em>et al.</em>, 2016).'
    return sentences.join(' ')
}

def toolBibliographyText(selected_methods) {
    def bibliography = [
        regenie: '<li>Mbatchou J, Barnard L, Backman J, et al. Computationally efficient whole-genome regression for quantitative and binary traits. <em>Nature Genetics</em>. 2021;53:1097-1103. doi: <a href="https://doi.org/10.1038/s41588-021-00870-7">10.1038/s41588-021-00870-7</a>.</li>',
        gcta_fastgwa: '<li>Jiang L, Zheng Z, Qi T, et al. A resource-efficient tool for mixed model association analysis of large-scale data. <em>Nature Genetics</em>. 2019;51:1749-1755. doi: <a href="https://doi.org/10.1038/s41588-019-0530-8">10.1038/s41588-019-0530-8</a>.</li>',
        gcta_greml: '<li>Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: a tool for genome-wide complex trait analysis. <em>American Journal of Human Genetics</em>. 2011;88:76-82. doi: <a href="https://doi.org/10.1016/j.ajhg.2010.11.011">10.1016/j.ajhg.2010.11.011</a>.</li>',
        gcta_greml_ldms: '<li>Yang J, Bakshi A, Zhu Z, et al. Genetic variance estimation with imputed variants finds negligible missing heritability for human height and body mass index. <em>Nature Genetics</em>. 2015;47:1114-1120. doi: <a href="https://doi.org/10.1038/ng.3390">10.1038/ng.3390</a>.</li>',
        gcta_bivariate_reml: '<li>Lee SH, Yang J, Goddard ME, Visscher PM, Wray NR. Estimation of pleiotropy between complex diseases using single-nucleotide polymorphism-derived genomic relationships and restricted maximum likelihood. <em>Bioinformatics</em>. 2012;28:2540-2542. doi: <a href="https://doi.org/10.1093/bioinformatics/bts474">10.1093/bioinformatics/bts474</a>.</li>',
        gcta_hereg: '<li>Chen GB. Estimating heritability of complex traits from genome-wide association studies using IBS-based Haseman-Elston regression. <em>Frontiers in Genetics</em>. 2014;5:107. doi: <a href="https://doi.org/10.3389/fgene.2014.00107">10.3389/fgene.2014.00107</a>.</li>',
        ldak_kvik: '<li>Hof JP, Speed D. LDAK-KVIK performs fast and powerful mixed-model association analysis of quantitative and binary phenotypes. <em>Nature Genetics</em>. 2025;57:2116-2123. doi: <a href="https://doi.org/10.1038/s41588-025-02286-z">10.1038/s41588-025-02286-z</a>.</li>',
        ldak: '<li>Speed D, Hemani G, Johnson MR, Balding DJ. Improved heritability estimation from genome-wide SNPs. <em>American Journal of Human Genetics</em>. 2012;91:1011-1021. doi: <a href="https://doi.org/10.1016/j.ajhg.2012.10.010">10.1016/j.ajhg.2012.10.010</a>.</li>',
        rhe_mc: '<li>Pazokitoroudi A, Wu Y, Burch KS, Hou K, Zhou A, Pasaniuc B, Sankararaman S. Efficient variance components analysis across millions of genomes. <em>Nature Communications</em>. 2020;11:4020. doi: <a href="https://doi.org/10.1038/s41467-020-17576-9">10.1038/s41467-020-17576-9</a>.</li>',
        mph: '<li>Jiang J. MPH: fast REML for large-scale genome partitioning of quantitative genetic variation. <em>Bioinformatics</em>. 2024;40:btae298. doi: <a href="https://doi.org/10.1093/bioinformatics/btae298">10.1093/bioinformatics/btae298</a>.</li>',
        ldak_sumstats: '<li>Speed D, Balding DJ. SumHer better estimates the SNP heritability of complex traits from summary statistics. <em>Nature Genetics</em>. 2019;51:277-284. doi: <a href="https://doi.org/10.1038/s41588-018-0279-5">10.1038/s41588-018-0279-5</a>.</li>',
        ldsc: '<li>Bulik-Sullivan BK, Loh PR, Finucane HK, et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. <em>Nature Genetics</em>. 2015;47:291-295. doi: <a href="https://doi.org/10.1038/ng.3211">10.1038/ng.3211</a>.</li>',
        gwaslab: '<li>GWASLab. <a href="https://cloufield.github.io/gwaslab/">https://cloufield.github.io/gwaslab/</a>.</li>',
        metasoft: '<li>Han B, Eskin E. Random-effects model aimed at discovering associations in meta-analysis of genome-wide association studies. <em>American Journal of Human Genetics</em>. 2011;88:586-598. doi: <a href="https://doi.org/10.1016/j.ajhg.2011.04.014">10.1016/j.ajhg.2011.04.014</a>.</li>',
        mrmega: '<li>Mägi R, Horikoshi M, Sofer T, et al. Trans-ethnic meta-regression of genome-wide association studies accounting for ancestry increases power for discovery and improves fine-mapping resolution. <em>Human Molecular Genetics</em>. 2017;26:3639-3650. doi: <a href="https://doi.org/10.1093/hmg/ddx280">10.1093/hmg/ddx280</a>.</li>',
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
