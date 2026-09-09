process LDAK_HE {
    tag "${meta.id}_${meta2.id}"
    label 'process_medium'
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
    tuple val(meta), path("${prefix}.he"), emit: he_results
    tuple val(meta), path("${prefix}.he.liab"), emit: he_liability, optional: true
    tuple val(meta), path("${prefix}.factor"), emit: factor, optional: true
    tuple val(meta), path("${prefix}.coeff"), emit: coeff, optional: true
    tuple val(meta), path("${prefix}.combined"), emit: combined, optional: true
    tuple val(meta), path("${prefix}.cross"), emit: cross, optional: true
    tuple val(meta), path("${prefix}.progress"), emit: progress, optional: true
    tuple val(meta), path("${prefix}.share"), emit: share, optional: true
    tuple val(meta), path("${prefix}.he.within"), emit: he_within, optional: true
    tuple val(meta), path("${prefix}.he.across"), emit: he_across, optional: true
    tuple val(meta), path("${prefix}.he.compare"), emit: he_compare, optional: true
    tuple val(meta), path("${prefix}.cross.within"), emit: cross_within, optional: true
    tuple val(meta), path("${prefix}.cross.across"), emit: cross_across, optional: true
    tuple val(meta), path("${prefix}.share.within"), emit: share_within, optional: true
    tuple val(meta), path("${prefix}.share.across"), emit: share_across, optional: true
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | grep -oP '(?<=^Version )[0-9.]+'"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prevalence_arg = prevalence ? "--prevalence \"${prevalence}\"" : ''
    def grm_prefix = grm_files.find { grm_file -> grm_file.name.endsWith('.grm.bin') }.name.replaceFirst(/\.grm\.bin$/, '')
    prefix = task.ext.prefix ?: "${meta.id}"
    def keep_arg = keep_file ? "--keep \"${keep_file}\"" : ''
    def quant_covar_arg = quant_covariates_file ? "--covar \"${quant_covariates_file}\"" : ''
    def cat_covar_arg = cat_covariates_file ? "--factors \"${cat_covariates_file}\"" : ''

    """
    ldak6 --he "${prefix}" \\
        --pheno "${phenotype_file}" \\
        --grm "${grm_prefix}" \\
        ${keep_arg} \\
        ${quant_covar_arg} \\
        ${cat_covar_arg} \\
        --max-threads "${task.cpus}" \\
        ${prevalence_arg} \\
        ${args} \\
        2>&1 | tee "${prefix}.log"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.he"
    touch "${prefix}.he.liab"
    touch "${prefix}.factor"
    touch "${prefix}.coeff"
    touch "${prefix}.combined"
    touch "${prefix}.cross"
    touch "${prefix}.progress"
    touch "${prefix}.share"
    touch "${prefix}.he.within"
    touch "${prefix}.he.across"
    touch "${prefix}.he.compare"
    touch "${prefix}.cross.within"
    touch "${prefix}.cross.across"
    touch "${prefix}.share.within"
    touch "${prefix}.share.across"
    touch "${prefix}.log"
    """
}
