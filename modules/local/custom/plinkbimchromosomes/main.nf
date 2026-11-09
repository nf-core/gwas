process CUSTOM_PLINKBIMCHROMOSOMES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(bim)

    output:
    tuple val(meta), path("${prefix}.chromosomes.txt"), emit: chromosomes
    tuple val("${task.process}"), val("coreutils"), eval("cat --version | head -n 1 | cut -d ' ' -f 4"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // `uniq` first collapses each contiguous run of one chromosome, so only a handful of lines reach the
    // order-preserving de-duplication that follows, however many variants the BIM holds.
    """
    sed 's/[[:space:]].*//' "${bim}" \\
        | uniq \\
        | cat -n \\
        | sort -s -u -k2,2 \\
        | sort -n -k1,1 \\
        | cut -f2 \\
        > "${prefix}.chromosomes.txt"
    """

    // The listing is the scatter's cardinality rather than a scientific result, and it costs one pass over
    // the BIM, so the stub computes it too: a stub run then scatters exactly as a real one would.
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    sed 's/[[:space:]].*//' "${bim}" \\
        | uniq \\
        | cat -n \\
        | sort -s -u -k2,2 \\
        | sort -n -k1,1 \\
        | cut -f2 \\
        > "${prefix}.chromosomes.txt"
    """
}
