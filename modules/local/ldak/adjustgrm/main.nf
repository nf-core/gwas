process LDAK_ADJUSTGRM {
    tag "${meta.id}_${meta2.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c94e5424f08cbd7e0856ab3c3a9992b4080944a5d0c497ce0abcceb413db4a3e/data'
        : 'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129'}"

    input:
    tuple val(meta), path(grm_files)
    tuple val(meta2), path(phenotype_file)
    tuple val(meta3), path(keep_file)
    tuple val(meta4), path(adjustment_covariates_file)

    output:
    tuple val(meta), path("${prefix}.grm.bin"), path("${prefix}.grm.id"), path("${prefix}.grm.details"), path("${prefix}.grm.adjust"), path("${prefix}.grm.root"), emit: adjusted_grm
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | grep -oP '(?<=^Version )[0-9.]+'"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def grm_prefix = grm_files.find { grm_file -> grm_file.name.endsWith('.grm.bin') }.name.replaceFirst(/\.grm\.bin$/, '')
    prefix = task.ext.prefix ?: "${meta.id}"
    def keep_arg = keep_file ? "--keep \"${keep_file}\"" : ''
    def covar_arg = adjustment_covariates_file ? "--covar \"${adjustment_covariates_file}\"" : ''
    """
    ldak6 --adjust-grm "${prefix}" \\
        --grm "${grm_prefix}" \\
        --pheno "${phenotype_file}" \\
        ${keep_arg} \\
        ${covar_arg} \\
        --max-threads "${task.cpus}" \\
        ${args} \\
        2>&1 | tee "${prefix}.log"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.grm.bin"
    touch "${prefix}.grm.id"
    touch "${prefix}.grm.details"
    touch "${prefix}.grm.adjust"
    touch "${prefix}.grm.root"
    touch "${prefix}.log"
    """
}
