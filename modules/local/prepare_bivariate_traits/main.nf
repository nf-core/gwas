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
    tuple val(meta), path("${prefix}.qcovar"), emit: quant_covariates, optional: true
    tuple val(meta), path("${prefix}.covar"), emit: cat_covariates, optional: true
    tuple val(meta), path("${prefix}.pair.log"), emit: log
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
    def qcovar_stub = pair_quant_covariates ? "printf 'stub stub 0\\n' > \"${prefix}.qcovar\"" : ''
    def covar_stub = pair_cat_covariates ? "printf 'stub stub 0\\n' > \"${prefix}.covar\"" : ''
    """
    printf 'stub stub 0 0\n' > "${prefix}.pheno"
    ${qcovar_stub}
    ${covar_stub}
    printf '%b\n' \
        'left_samples\t1' \
        'right_samples\t1' \
        'endpoint_overlap_samples\t1' \
        'union_samples\t1' \
        'left_nonmissing\t1' \
        'right_nonmissing\t1' \
        'both_nonmissing\t1' \
        'quantitative_covariate_samples\t${pair_quant_covariates ? 1 : 0}' \
        'categorical_covariate_samples\t${pair_cat_covariates ? 1 : 0}' \
        > "${prefix}.pair.log"
    printf '"%s":\n    python: %s\n' \
        '${task.process}' \
        "\$(python3 --version | sed 's/^Python //')" \
        > versions.yml
    """
}
