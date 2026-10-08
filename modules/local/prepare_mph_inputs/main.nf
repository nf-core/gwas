process PREPARE_MPH_INPUTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    tuple val(meta), path(phenotype_table), path(quant_covariates), path(cat_covariates), val(trait_names), val(missing_tokens)
    tuple val(meta2), path(grm_iid), path(fam)

    output:
    tuple val(meta), path("${prefix}.mph.pheno.csv"), emit: phenotype
    tuple val(meta), path("${prefix}.mph.covar.csv"), emit: covariates, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // `prefix` must remain visible to the output declarations. Every other assignment must remain visible to
    // the template and hold a JSON string literal, which is also valid Python syntax.
    prefix = task.ext.prefix ?: "${meta.id}"
    phenotype_table_literal = groovy.json.JsonOutput.toJson(phenotype_table.toString())
    quant_covariates_literal = groovy.json.JsonOutput.toJson(quant_covariates ? quant_covariates.toString() : '')
    cat_covariates_literal = groovy.json.JsonOutput.toJson(cat_covariates ? cat_covariates.toString() : '')
    grm_iid_literal = groovy.json.JsonOutput.toJson(grm_iid.toString())
    fam_literal = groovy.json.JsonOutput.toJson(fam.toString())
    prefix_literal = groovy.json.JsonOutput.toJson(prefix.toString())
    analysis_id_literal = groovy.json.JsonOutput.toJson(meta.id.toString())
    task_process_literal = groovy.json.JsonOutput.toJson(task.process.toString())
    // Serialised twice on purpose: the inner call renders the structure as JSON, the outer one renders that
    // JSON as a quoted string. A bare JSON object is not valid Python -- its `null`, `true` and `false` are
    // not Python literals -- so the template parses a string rather than embedding an expression.
    trait_names_literal = groovy.json.JsonOutput.toJson(groovy.json.JsonOutput.toJson(trait_names))
    missing_tokens_literal = groovy.json.JsonOutput.toJson(groovy.json.JsonOutput.toJson(missing_tokens))
    template('prepare_mph_inputs.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def stub_traits = trait_names.join(',')
    def stub_trait_values = trait_names.collect { '0' }.join(',')
    // The covariate file is written only when the analysis supplied covariates, so a stub run exercises the
    // same optional-output cardinality the real script produces.
    def stub_covariates = quant_covariates || cat_covariates
        ? "printf 'IID,intercept,stub_covariate\\nstub,1,0\\n' > \"${prefix}.mph.covar.csv\""
        : ''
    """
    printf 'IID,${stub_traits}\\nstub,${stub_trait_values}\\n' > "${prefix}.mph.pheno.csv"
    ${stub_covariates}
    printf '"%s":\\n    python: %s\\n' \\
        '${task.process}' \\
        "\$(python3 --version | sed 's/^Python //')" \\
        > versions.yml
    """
}
