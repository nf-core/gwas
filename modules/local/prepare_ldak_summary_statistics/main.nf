process PREPARE_LDAK_SUMMARY_STATISTICS {
    tag "${meta.summary_statistics_id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    tuple val(meta), path(gwaslab_summary_statistics)

    output:
    tuple val(meta), path("${prefix}.summaries"), emit: summary_statistics
    path 'versions.yml', emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.summary_statistics_id}.ldak"
    input_literal = groovy.json.JsonOutput.toJson(gwaslab_summary_statistics.toString())
    output_literal = groovy.json.JsonOutput.toJson("${prefix}.summaries")
    task_process_literal = groovy.json.JsonOutput.toJson(task.process.toString())
    template('prepare_ldak_summary_statistics.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.summary_statistics_id}.ldak"
    """
    printf 'Predictor\tA1\tA2\tZ\tn\tA1Freq\n' > "${prefix}.summaries"
    printf '"${task.process}":\n    python: %s\n' "\$(python3 --version | sed 's/^Python //')" > versions.yml
    """
}
