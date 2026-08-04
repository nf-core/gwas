process NORMALISE_PHENOTYPES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    // Staged under `input/` rather than at the task root, which is not cosmetic. Nextflow stages an
    // input as a symlink into the task directory and places no guard on an output whose name equals
    // a staged input's, so an analysis unit named after its own phenotype file — `cohort1` beside
    // `cohort1.pheno`, which is exactly how the fixture bundle names things — would have this module
    // write `cohort1.pheno` straight through the symlink and destroy the researcher's source file,
    // silently, with the run reporting success. This module is the only one in the repo that reads
    // and writes the same extensions, so it is the only one that can collide; `input/` makes the
    // collision impossible rather than merely unlikely.
    tuple val(meta), path(phenotype, stageAs: 'input/*'), path(quant_covariates, stageAs: 'input/*'), path(cat_covariates, stageAs: 'input/*')

    output:
    tuple val(meta), path("${prefix}.pheno"), emit: phenotype
    tuple val(meta), path("${prefix}.noheader.pheno"), emit: phenotype_headerless
    tuple val(meta), path("${prefix}.qcovar"), emit: quant_covariates, optional: true
    tuple val(meta), path("${prefix}.noheader.qcovar"), emit: quant_covariates_headerless, optional: true
    tuple val(meta), path("${prefix}.catcovar"), emit: cat_covariates, optional: true
    tuple val(meta), path("${prefix}.noheader.catcovar"), emit: cat_covariates_headerless, optional: true
    tuple val(meta), path("${prefix}.covar"), emit: covariates, optional: true
    tuple val(meta), path("${prefix}.adjustcovar"), emit: adjustment_covariates, optional: true
    tuple val(meta), path("${prefix}.normalise.log"), emit: log
    // `eval()` is the house style for version capture, but Nextflow rejects an `eval` output on any
    // process whose script is not Bash — and this one is a Python `template`. So the interpreter
    // version is written from inside the template instead, exactly as the repo's other template
    // module `modules/local/ldak/calcinflation` does.
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // `prefix` must remain visible to the output declarations. The other assignments must remain
    // visible to the template and contain JSON string literals, which are also valid Python syntax.
    prefix = task.ext.prefix ?: "${meta.id}"
    phenotype_literal = groovy.json.JsonOutput.toJson(phenotype.toString())
    quant_covariates_literal = groovy.json.JsonOutput.toJson(quant_covariates ? quant_covariates.toString() : '')
    cat_covariates_literal = groovy.json.JsonOutput.toJson(cat_covariates ? cat_covariates.toString() : '')
    phenotype_column_literal = groovy.json.JsonOutput.toJson(meta.phenotype_column.toString())
    trait_type_literal = groovy.json.JsonOutput.toJson(meta.is_binary ? 'binary' : 'quantitative')
    case_value_literal = groovy.json.JsonOutput.toJson((meta.case_value ?: '').toString())
    control_value_literal = groovy.json.JsonOutput.toJson((meta.control_value ?: '').toString())
    prefix_literal = groovy.json.JsonOutput.toJson(prefix.toString())
    analysis_id_literal = groovy.json.JsonOutput.toJson(meta.id.toString())
    task_process_literal = groovy.json.JsonOutput.toJson(task.process.toString())
    template('normalise_phenotypes.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def quant_covariates_stub = quant_covariates
        ? """
    touch "${prefix}.qcovar"
    touch "${prefix}.noheader.qcovar"
    """
        : ''
    def cat_covariates_stub = cat_covariates
        ? """
    touch "${prefix}.catcovar"
    touch "${prefix}.noheader.catcovar"
    """
        : ''
    def merged_covariates_stub = quant_covariates || cat_covariates ? """touch "${prefix}.covar"\n""" : ''
    def adjustment_covariates_stub = quant_covariates || cat_covariates ? """touch "${prefix}.adjustcovar"\n""" : ''
    """
    printf 'FID\\tIID\\tPHENO\\n' > "${prefix}.pheno"
    printf '' > "${prefix}.noheader.pheno"
    ${quant_covariates_stub}
    ${cat_covariates_stub}
    ${merged_covariates_stub}
    ${adjustment_covariates_stub}
    touch "${prefix}.normalise.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """
}
