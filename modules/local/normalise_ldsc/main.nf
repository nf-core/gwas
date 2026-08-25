process NORMALISE_LDSC {
    tag "${meta.request_id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    tuple val(meta), path(observed_log), path(liability_log), path(munging_logs)

    output:
    tuple val(meta), path('heritability.tsv'), emit: heritability
    tuple val(meta), path('genetic_correlation.tsv'), optional: true, emit: genetic_correlation
    tuple val(meta), path('genetic_covariance.tsv'), optional: true, emit: genetic_covariance
    tuple val(meta), path('diagnostics.tsv'), emit: diagnostics
    tuple val(meta), path('provenance.json'), emit: provenance
    path 'versions.yml', emit: versions, topic: versions

    script:
    meta_literal = groovy.json.JsonOutput.toJson(groovy.json.JsonOutput.toJson(meta))
    observed_log_literal = groovy.json.JsonOutput.toJson(observed_log.toString())
    liability_log_literal = liability_log ? groovy.json.JsonOutput.toJson(liability_log.toString()) : 'None'
    munging_logs_literal = groovy.json.JsonOutput.toJson(munging_logs.collect { it.toString() })
    task_process_literal = groovy.json.JsonOutput.toJson(task.process.toString())
    task_container_literal = groovy.json.JsonOutput.toJson(task.container.toString())
    ldsc_container_literal = groovy.json.JsonOutput.toJson(task.ext.ldsc_container.toString())
    template('normalise_ldsc.py')

    stub:
    def pairwise = meta.method == 'ldsc_rg'
    def endpoint = pairwise ? 'left' : 'unary'
    def summary_id = pairwise ? meta.left_summary_statistics_id : meta.summary_statistics_id
    def trait_id = pairwise ? meta.left_trait_id : meta.trait_id
    def trait_type = pairwise ? meta.left_trait_type : meta.trait_type
    """
    cat <<-'END_HERITABILITY' > heritability.tsv
    relationship_id\trequest_id\tmethod\tendpoint\tsummary_statistics_id\ttrait_id\ttrait_type\tscale\testimate\tstandard_error\tclassification\tnative_artifact
    ${meta.relationship_id ?: 'NA'}\t${meta.request_id}\t${meta.method}\t${endpoint}\t${summary_id}\t${trait_id}\t${trait_type}\tobserved\tNA\tNA\tcompleted_nonestimable\trequests/${meta.method}/${meta.request_id}/native.observed.log
    END_HERITABILITY
    ${pairwise ? """
    cat <<-'END_CORRELATION' > genetic_correlation.tsv
    relationship_id\trequest_id\tmethod\tleft_summary_statistics_id\tright_summary_statistics_id\tleft_trait_id\tright_trait_id\testimate\tstandard_error\tz_score\tp_value\tclassification\tnative_artifact
    ${meta.relationship_id}\t${meta.request_id}\t${meta.method}\t${meta.left_summary_statistics_id}\t${meta.right_summary_statistics_id}\t${meta.left_trait_id}\t${meta.right_trait_id}\tNA\tNA\tNA\tNA\tcompleted_nonestimable\trequests/${meta.method}/${meta.request_id}/native.observed.log
    END_CORRELATION
    cat <<-'END_COVARIANCE' > genetic_covariance.tsv
    relationship_id\trequest_id\tmethod\tleft_summary_statistics_id\tright_summary_statistics_id\tleft_trait_id\tright_trait_id\tscale\testimate\tstandard_error\tclassification\tnative_artifact
    ${meta.relationship_id}\t${meta.request_id}\t${meta.method}\t${meta.left_summary_statistics_id}\t${meta.right_summary_statistics_id}\t${meta.left_trait_id}\t${meta.right_trait_id}\tobserved\tNA\tNA\tcompleted_nonestimable\trequests/${meta.method}/${meta.request_id}/native.observed.log
    END_COVARIANCE
    """ : ''}
    cat <<-'END_DIAGNOSTICS' > diagnostics.tsv
    relationship_id\trequest_id\tmetric\tvalue
    ${meta.relationship_id ?: 'NA'}\t${meta.request_id}\tclassification\tcompleted_nonestimable
    END_DIAGNOSTICS
    cat <<-'END_PROVENANCE' > provenance.json
    {"schema_version":"1.0","request_id":"${meta.request_id}","method":"${meta.method}","classification":"completed_nonestimable","stub":true}
    END_PROVENANCE
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """
}
