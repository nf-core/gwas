process LDAK_FASTPCGC {
    tag "${meta.id}_${meta2.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c94e5424f08cbd7e0856ab3c3a9992b4080944a5d0c497ce0abcceb413db4a3e/data'
        : 'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129'}"

    input:
    tuple val(meta), path(phenotype_file), val(prevalence)
    tuple val(meta2), path(bed), path(bim), path(fam), val(power)
    tuple val(meta3), path(weights_file)
    tuple val(meta4), path(quant_covariates_file)
    tuple val(meta5), path(cat_covariates_file)

    output:
    tuple val(meta), path("${prefix}.fastpcgc"), emit: fastpcgc_results
    tuple val(meta), path("${prefix}.fastpcgc.marginal"), emit: fastpcgc_marginal, optional: true
    tuple val(meta), path("${prefix}.cats"), emit: cats, optional: true
    tuple val(meta), path("${prefix}.coeff"), emit: coeff, optional: true
    tuple val(meta), path("${prefix}.combined"), emit: combined, optional: true
    tuple val(meta), path("${prefix}.cross"), emit: cross, optional: true
    tuple val(meta), path("${prefix}.enrich"), emit: enrich, optional: true
    tuple val(meta), path("${prefix}.share"), emit: share, optional: true
    tuple val(meta), path("${prefix}.ind.hers"), emit: ind_hers, optional: true
    tuple val(meta), path("${prefix}.jackests"), emit: jackests, optional: true
    tuple val(meta), path("${prefix}.labels"), emit: labels, optional: true
    tuple val(meta), path("${prefix}.progress"), emit: progress, optional: true
    tuple val(meta), path("${prefix}.repetitions"), emit: repetitions, optional: true
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | grep -oP '(?<=^Version )[0-9.]+'"), emit: versions_ldak6, topic: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def weights_arg = weights_file ? "--weights \"${weights_file}\"" : ''
    def quant_covar_arg = quant_covariates_file ? "--covar \"${quant_covariates_file}\"" : ''
    def cat_covar_arg = cat_covariates_file ? "--factors \"${cat_covariates_file}\"" : ''
    """
    ldak6 --fast-pcgc "${prefix}" \\
        --bfile "${bed.baseName}" \\
        --pheno "${phenotype_file}" \\
        --power "${power}" \\
        ${weights_arg} \\
        ${quant_covar_arg} \\
        ${cat_covar_arg} \\
        --prevalence "${prevalence}" \\
        --max-threads "${task.cpus}" \\
        ${args} \\
        2>&1 | tee "${prefix}.log"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.fastpcgc"
    touch "${prefix}.fastpcgc.marginal"
    touch "${prefix}.cats"
    touch "${prefix}.coeff"
    touch "${prefix}.combined"
    touch "${prefix}.cross"
    touch "${prefix}.enrich"
    touch "${prefix}.share"
    touch "${prefix}.ind.hers"
    touch "${prefix}.jackests"
    touch "${prefix}.labels"
    touch "${prefix}.progress"
    touch "${prefix}.repetitions"
    touch "${prefix}.log"
    """
}
