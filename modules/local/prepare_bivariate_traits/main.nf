process PREPARE_BIVARIATE_TRAITS {
    tag "${meta.request_id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    tuple val(meta), path(left_phenotype, stageAs: 'left/*'), path(right_phenotype, stageAs: 'right/*'), path(pair_quant_covariates, stageAs: 'pair_qcov/*'), path(pair_cat_covariates, stageAs: 'pair_cov/*')

    output:
    tuple val(meta), path("${prefix}.pheno"), emit: phenotype
    // One validated pair covariate table, serialised twice for the two interfaces that read it. GCTA rejects a
    // header row, while MPH names its covariates on the command line and reports them back by name, so the
    // headerless files feed the GCTA estimators and the headered ones feed the MPH serializer. Both are written
    // from the same validated rows, so the researcher's file is normalised exactly once and no estimator
    // independently redefines the pair's covariate set. The suffixes are `prepare_phenotype_inputs`': headered
    // is `.qcovar`/`.catcovar` and headerless is `.noheader.*`.
    tuple val(meta), path("${prefix}.qcovar"), emit: named_quant_covariates, optional: true
    tuple val(meta), path("${prefix}.catcovar"), emit: named_cat_covariates, optional: true
    tuple val(meta), path("${prefix}.noheader.qcovar"), emit: quant_covariates, optional: true
    tuple val(meta), path("${prefix}.noheader.catcovar"), emit: cat_covariates, optional: true
    path 'versions.yml', emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.request_id
    left_phenotype_literal = groovy.json.JsonOutput.toJson(left_phenotype.toString())
    right_phenotype_literal = groovy.json.JsonOutput.toJson(right_phenotype.toString())
    pair_quant_covariates_literal = groovy.json.JsonOutput.toJson(pair_quant_covariates ? pair_quant_covariates.toString() : '')
    pair_cat_covariates_literal = groovy.json.JsonOutput.toJson(pair_cat_covariates ? pair_cat_covariates.toString() : '')
    prefix_literal = groovy.json.JsonOutput.toJson(prefix.toString())
    request_id_literal = groovy.json.JsonOutput.toJson(meta.request_id.toString())
    task_process_literal = groovy.json.JsonOutput.toJson(task.process.toString())
    template('prepare_bivariate_traits.py')

    stub:
    prefix = task.ext.prefix ?: meta.request_id
    def qcovar_stub = pair_quant_covariates
        ? "printf 'FID IID STUBQ\\n' > \"${prefix}.qcovar\"\n    printf 'stub stub 0\\n' > \"${prefix}.noheader.qcovar\""
        : ''
    def covar_stub = pair_cat_covariates
        ? "printf 'FID IID STUBC\\n' > \"${prefix}.catcovar\"\n    printf 'stub stub 0\\n' > \"${prefix}.noheader.catcovar\""
        : ''
    """
    printf 'stub stub 0 0\n' > "${prefix}.pheno"
    ${qcovar_stub}
    ${covar_stub}
    printf '"%s":\n    python: %s\n' \
        '${task.process}' \
        "\$(python3 --version | sed 's/^Python //')" \
        > versions.yml
    """
}
