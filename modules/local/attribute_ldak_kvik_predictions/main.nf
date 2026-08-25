process ATTRIBUTE_LDAK_KVIK_PREDICTIONS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(root), path(loco_details), path(loco_prs)

    output:
    tuple val(meta), path("${prefix}.ldak_kvik.step1.root"), path("${prefix}.ldak_kvik.step1.loco.details"), path("${prefix}.ldak_kvik.step1.loco.prs"), emit: bundle
    tuple val("${task.process}"), val("coreutils"), eval("cp --version | head -n 1 | cut -d ' ' -f 4"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp "${root}" "${prefix}.ldak_kvik.step1.root"
    cp "${loco_details}" "${prefix}.ldak_kvik.step1.loco.details"
    cp "${loco_prs}" "${prefix}.ldak_kvik.step1.loco.prs"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.ldak_kvik.step1.root"
    touch "${prefix}.ldak_kvik.step1.loco.details"
    touch "${prefix}.ldak_kvik.step1.loco.prs"
    """
}
