process CUSTOM_GCTACONCATLDSCORES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(ld_scores, arity: '1..*')

    output:
    tuple val(meta), path("${prefix}_gcta_ld.score.ld"), emit: ld_scores
    tuple val("${task.process}"), val("coreutils"), eval("cat --version | head -n 1 | cut -d ' ' -f 4"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def tables = ld_scores.collect { table -> "\"${table}\"" }.join(' ')
    """
    head -n 1 "${ld_scores[0]}" > "${prefix}_gcta_ld.score.ld"
    for table in ${tables}; do
        tail -n +2 "\${table}" >> "${prefix}_gcta_ld.score.ld"
    done
    """

    // Concatenation costs no more than copying the stub tables it is given, so the stub runs it unchanged.
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def tables = ld_scores.collect { table -> "\"${table}\"" }.join(' ')
    """
    head -n 1 "${ld_scores[0]}" > "${prefix}_gcta_ld.score.ld"
    for table in ${tables}; do
        tail -n +2 "\${table}" >> "${prefix}_gcta_ld.score.ld"
    done
    """
}
