process GCTA_BIVARIATEHEREGLDMS {
    tag "${meta.id}_${meta2.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46b0d05f0daa47561d87d2a9cac5e51edc2c78e26f1bbab439c688386241a274/data'
        : 'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9'}"

    input:
    tuple val(meta), path(grm_files), val(grm_prefixes)
    tuple val(meta2), path(phenotypes_file), val(phenotype_col1), val(phenotype_col2)

    output:
    tuple val(meta), path("*.HEreg"), emit: hereg_results
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val("gcta"), eval("gcta --version | sed -En 's/^[*] version v([0-9.]*).*/\\1/p'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mgrm_entries = grm_prefixes.collect { grm_prefix -> "\"${grm_prefix}\"" }.join(' ')
    def mgrm_file = "${prefix}.mgrm"
    def hereg_bivar_param = phenotype_col1 && phenotype_col2 ? "--HEreg-bivar ${phenotype_col1} ${phenotype_col2}" : '--HEreg-bivar'
    """
    printf '%s\n' ${mgrm_entries} > "${mgrm_file}"

    gcta \\
        ${hereg_bivar_param} \\
        --mgrm "${mgrm_file}" \\
        --pheno "${phenotypes_file}" \\
        --out "${prefix}" \\
        --thread-num "${task.cpus}" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.HEreg"
    touch "${prefix}.log"
    """
}
