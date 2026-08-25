process NORMALISE_LDAK_SUMHER {
    tag "${meta.request_id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    tuple val(meta), path(hers), path(extra), path(overlap), path(ldak_log), path(preparation), path(hers_liability)

    output:
    tuple val(meta), path('heritability.tsv'), emit: heritability
    tuple val(meta), path('diagnostics.tsv'), emit: diagnostics
    tuple val(meta), path('provenance.json'), emit: provenance
    path 'versions.yml', emit: versions, topic: versions

    script:
    meta_literal = groovy.json.JsonOutput.toJson(groovy.json.JsonOutput.toJson(meta))
    hers_literal = groovy.json.JsonOutput.toJson(hers.toString())
    extra_literal = groovy.json.JsonOutput.toJson(extra.toString())
    overlap_literal = groovy.json.JsonOutput.toJson(overlap.toString())
    log_literal = groovy.json.JsonOutput.toJson(ldak_log.toString())
    preparation_literal = groovy.json.JsonOutput.toJson(preparation.toString())
    hers_liability_literal = groovy.json.JsonOutput.toJson(hers_liability ? hers_liability.toString() : '')
    task_process_literal = groovy.json.JsonOutput.toJson(task.process.toString())
    task_container_literal = groovy.json.JsonOutput.toJson(task.container.toString())
    template('normalise_ldak_sumher.py')

    stub:
    """
    printf '%s\n' \
        'summary_statistics_id\trequest_id\tmethod\ttrait_id\ttrait_type\tscale\testimate\tstandard_error\tclassification' \
        '${meta.summary_statistics_id}\t${meta.request_id}\t${meta.method}\t${meta.trait_id}\t${meta.trait_type}\tobserved\tNA\tNA\tcompleted_nonestimable' \
        > heritability.tsv
    printf '%s\n' \
        'request_id\tmetric\tvalue' \
        '${meta.request_id}\tclassification\tcompleted_nonestimable' \
        > diagnostics.tsv
    printf '{"schema_version":"1.0","request_id":"%s","method":"%s","classification":"completed_nonestimable","stub":true}\n' '${meta.request_id}' '${meta.method}' > provenance.json
    printf '"${task.process}":\n    python: %s\n' "\$(python3 --version | sed 's/^Python //')" > versions.yml
    """
}
