process PREPARE_GENOTYPE_SOURCE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    // The supplied genotype bundle, staged under `source/` so no output of this process can ever be written
    // through a staged symlink onto the researcher's own file. The members arrive in manifest-group order and
    // the emitted table preserves that order, so the caller supplies the order and this module never repeats
    // the manifest's genotype-group contract.
    tuple val(meta), path(genotype_files, stageAs: 'source/*')

    output:
    // One row per member: role, basename, byte size, SHA-256 of the bytes. This is provenance, and it is
    // deliberately a separate output from the identity below so that a caller can keep member names out of a
    // task hash while still publishing them.
    tuple val(meta), path("${prefix}.genotype_source.tsv"), emit: source_members
    // The bundle's content identity: one SHA-256 over the role/digest columns, so two bundles whose members
    // hold the same bytes under the same roles are identical whatever the researcher named the files.
    tuple val(meta), path("${prefix}.genotype_source.sha256"), emit: source_identity
    // An `eval` output is only accepted for a Bash script, and this process is a Python template, so the
    // version row is written by the template itself, as every other template module here does.
    path 'versions.yml', emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // `prefix` must remain visible to the output declarations above.
    prefix = task.ext.prefix ?: "${meta.id}"
    // The double encoding renders a Python `str` the template decodes with `json.loads`, matching the
    // convention every other template module in this repository uses.
    members_literal = groovy.json.JsonOutput.toJson(groovy.json.JsonOutput.toJson((genotype_files instanceof List ? genotype_files : [genotype_files]).collect { member -> member.toString() }))
    prefix_literal = groovy.json.JsonOutput.toJson(prefix.toString())
    task_process_literal = groovy.json.JsonOutput.toJson(task.process.toString())
    template('prepare_genotype_source.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Both outputs are unconditional, so the stub emits the same two files for every source format. A stub
    // whose cardinality varied by format would make a stubbed subworkflow test unable to see one shape.
    """
    printf 'bed\\tstub.bed\\t0\\t%064d\\n' 0 > "${prefix}.genotype_source.tsv"
    printf '%064d\\n' 0 > "${prefix}.genotype_source.sha256"
    printf '"%s":\\n    python: %s\\n' \\
        '${task.process}' \\
        "\$(python3 --version | sed 's/^Python //')" \\
        > versions.yml
    """
}
