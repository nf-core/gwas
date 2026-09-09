process CUSTOM_GCTAMERGEGRMPARTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(grm_files), val(ordered_parts)

    output:
    tuple val(meta), path("${prefix}.grm.*"), emit: grm_files
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def grm_bin_parts = ordered_parts.collect { part -> "\"${part.bin}\"" }.join(' ')
    def grm_n_bin_parts = ordered_parts.collect { part -> "\"${part.n_bin}\"" }.join(' ')
    def grm_id_parts = ordered_parts.collect { part -> "\"${part.id}\"" }.join(' ')
    """
    cat ${grm_bin_parts} > "${prefix}.grm.bin"
    cat ${grm_n_bin_parts} > "${prefix}.grm.N.bin"
    cat ${grm_id_parts} > "${prefix}.grm.id"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(sort --version | sed -n '1{s/sort (GNU coreutils) //;p}')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.grm.bin"
    touch "${prefix}.grm.N.bin"
    touch "${prefix}.grm.id"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(cat --version | head -n 1 | cut -d ' ' -f 4)
    END_VERSIONS
    """
}
