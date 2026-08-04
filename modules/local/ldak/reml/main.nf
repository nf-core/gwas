process LDAK_REML {
    tag "${meta.id}"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c94e5424f08cbd7e0856ab3c3a9992b4080944a5d0c497ce0abcceb413db4a3e/data'
        : 'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129'}"

    input:
    tuple val(meta), path(phenotype_file), val(prevalence)
    tuple val(meta2), path(grm_files)
    tuple val(meta3), path(keep_file)
    tuple val(meta4), path(quant_covariates_file)
    tuple val(meta5), path(cat_covariates_file)

    output:
    tuple val(meta), path("*.reml"), emit: reml_results
    tuple val(meta), path("*.coeff"), emit: coeff
    tuple val(meta), path("*.combined"), emit: combined, optional: true
    tuple val(meta), path("*.cross"), emit: cross
    tuple val(meta), path("*.indi.blp"), emit: indi_blp
    tuple val(meta), path("*.indi.res"), emit: indi_res
    tuple val(meta), path("*.progress"), emit: progress
    tuple val(meta), path("*.share"), emit: share
    tuple val(meta), path("*.vars"), emit: vars
    tuple val(meta), path("*.reml.liab"), emit: reml_liability, optional: true
    tuple val(meta), path("*.coeff.liab"), emit: coeff_liability, optional: true
    tuple val(meta), path("*.indi.blp.liab"), emit: indi_blp_liability, optional: true
    tuple val(meta), path("*.factor"), emit: liability_factor, optional: true
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | grep -oP '(?<=^Version )[0-9.]+'"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def grm_prefix = grm_files.find { grm_file -> grm_file.name.endsWith('.grm.bin') }.name.replaceFirst(/\.grm\.bin$/, '')
    def prefix = task.ext.prefix ?: meta.id
    def keep_arg = keep_file ? "--keep \"${keep_file}\"" : ''
    def quant_covar_arg = quant_covariates_file ? "--covar \"${quant_covariates_file}\"" : ''
    def cat_covar_arg = cat_covariates_file ? "--factors \"${cat_covariates_file}\"" : ''
    def prevalence_arg = prevalence ? "--prevalence \"${prevalence}\"" : ''

    """
    ldak6 --reml "${prefix}" \\
        --pheno "${phenotype_file}" \\
        --grm "${grm_prefix}" \\
        ${keep_arg} \\
        ${quant_covar_arg} \\
        ${cat_covar_arg} \\
        ${prevalence_arg} \\
        --max-threads "${task.cpus}" \\
        ${args} \\
        2>&1 | tee "${prefix}.log"
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    def liability_outputs = prevalence
        ? """
    touch "${prefix}.reml.liab"
    touch "${prefix}.coeff.liab"
    touch "${prefix}.indi.blp.liab"
    touch "${prefix}.factor"
    """
        : ''
    """
    touch "${prefix}.reml"
    touch "${prefix}.coeff"
    touch "${prefix}.combined"
    touch "${prefix}.cross"
    touch "${prefix}.indi.blp"
    touch "${prefix}.indi.res"
    touch "${prefix}.progress"
    touch "${prefix}.share"
    touch "${prefix}.vars"
    touch "${prefix}.log"
    ${liability_outputs}
    """
}
