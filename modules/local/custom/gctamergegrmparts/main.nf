process CUSTOM_GCTAMERGEGRMPARTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(grm_bin_parts), path(grm_n_bin_parts), path(grm_id_parts)

    output:
    tuple val(meta), path("${prefix}.grm.*"), emit: grm_files
    tuple val("${task.process}"), val("coreutils"), eval("cat --version | head -n 1 | cut -d ' ' -f 4"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def bin_parts = grm_bin_parts.collect { part_file -> "\"${part_file}\"" }.join(' ')
    def n_bin_parts = grm_n_bin_parts.collect { part_file -> "\"${part_file}\"" }.join(' ')
    def id_parts = grm_id_parts.collect { part_file -> "\"${part_file}\"" }.join(' ')
    """
    cat ${bin_parts} > "${prefix}.grm.bin"
    cat ${n_bin_parts} > "${prefix}.grm.N.bin"
    cat ${id_parts} > "${prefix}.grm.id"

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.grm.bin"
    touch "${prefix}.grm.N.bin"
    touch "${prefix}.grm.id"

    """
}
