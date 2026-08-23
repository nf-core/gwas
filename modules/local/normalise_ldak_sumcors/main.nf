process NORMALISE_LDAK_SUMCORS {
    tag "${meta.request_id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    tuple val(meta), path(correlations), path(correlations_full), path(overlap), path(ldak_log), path(left_preparation), path(right_preparation), path(correlations_liability)

    output:
    tuple val(meta), path('heritability.tsv'), emit: heritability
    tuple val(meta), path('genetic_correlation.tsv'), emit: genetic_correlation
    tuple val(meta), path('genetic_covariance.tsv'), emit: genetic_covariance
    tuple val(meta), path('diagnostics.tsv'), emit: diagnostics
    tuple val(meta), path('provenance.json'), emit: provenance
    path 'versions.yml', emit: versions, topic: versions

    script:
    meta_literal = groovy.json.JsonOutput.toJson(groovy.json.JsonOutput.toJson(meta))
    correlations_literal = groovy.json.JsonOutput.toJson(correlations.toString())
    correlations_full_literal = groovy.json.JsonOutput.toJson(correlations_full.toString())
    overlap_literal = groovy.json.JsonOutput.toJson(overlap.toString())
    log_literal = groovy.json.JsonOutput.toJson(ldak_log.toString())
    left_preparation_literal = groovy.json.JsonOutput.toJson(left_preparation.toString())
    right_preparation_literal = groovy.json.JsonOutput.toJson(right_preparation.toString())
    correlations_liability_literal = groovy.json.JsonOutput.toJson(correlations_liability ? correlations_liability.toString() : '')
    task_process_literal = groovy.json.JsonOutput.toJson(task.process.toString())
    task_container_literal = groovy.json.JsonOutput.toJson(task.container.toString())
    template('normalise_ldak_sumcors.py')

    stub:
    """
    printf '%s\n' \
        'relationship_id\trequest_id\tmethod\tendpoint\tsummary_statistics_id\ttrait_id\ttrait_type\tscale\testimate\tstandard_error\tclassification' \
        '${meta.relationship_id}\t${meta.request_id}\t${meta.method}\tleft\t${meta.left_summary_statistics_id}\t${meta.left_trait_id}\t${meta.left_trait_type}\tobserved\tNA\tNA\tcompleted_nonestimable' \
        '${meta.relationship_id}\t${meta.request_id}\t${meta.method}\tright\t${meta.right_summary_statistics_id}\t${meta.right_trait_id}\t${meta.right_trait_type}\tobserved\tNA\tNA\tcompleted_nonestimable' \
        > heritability.tsv
    printf '%s\n' \
        'relationship_id\trequest_id\tmethod\tleft_summary_statistics_id\tright_summary_statistics_id\tleft_trait_id\tright_trait_id\testimate\tstandard_error\tclassification' \
        '${meta.relationship_id}\t${meta.request_id}\t${meta.method}\t${meta.left_summary_statistics_id}\t${meta.right_summary_statistics_id}\t${meta.left_trait_id}\t${meta.right_trait_id}\tNA\tNA\tcompleted_nonestimable' \
        > genetic_correlation.tsv
    printf '%s\n' \
        'relationship_id\trequest_id\tmethod\tleft_summary_statistics_id\tright_summary_statistics_id\tleft_trait_id\tright_trait_id\tscale\testimate\tstandard_error\tclassification' \
        '${meta.relationship_id}\t${meta.request_id}\t${meta.method}\t${meta.left_summary_statistics_id}\t${meta.right_summary_statistics_id}\t${meta.left_trait_id}\t${meta.right_trait_id}\tobserved\tNA\tNA\tcompleted_nonestimable' \
        > genetic_covariance.tsv
    printf '%s\n' \
        'request_id\tmetric\tvalue' \
        '${meta.request_id}\tclassification\tcompleted_nonestimable' \
        > diagnostics.tsv
    printf '{"schema_version":"1.0","relationship_id":"%s","request_id":"%s","method":"%s","classification":"completed_nonestimable","stub":true}\n' '${meta.relationship_id}' '${meta.request_id}' '${meta.method}' > provenance.json
    printf '"${task.process}":\n    python: %s\n' "\$(python3 --version | sed 's/^Python //')" > versions.yml
    """
}
