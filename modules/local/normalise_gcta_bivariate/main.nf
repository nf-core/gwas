process NORMALISE_GCTA_BIVARIATE {
    tag "${meta.request_id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    tuple val(meta), path(hsq), path(gcta_log), path(pair_log)

    output:
    tuple val(meta), path('heritability.tsv'), emit: heritability
    tuple val(meta), path('genetic_correlation.tsv'), emit: genetic_correlation
    tuple val(meta), path('genetic_covariance.tsv'), emit: genetic_covariance
    tuple val(meta), path('diagnostics.tsv'), emit: diagnostics
    tuple val(meta), path('provenance.json'), emit: provenance
    path 'versions.yml', emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    meta_literal = groovy.json.JsonOutput.toJson(groovy.json.JsonOutput.toJson(meta))
    hsq_literal = groovy.json.JsonOutput.toJson(hsq.toString())
    gcta_log_literal = groovy.json.JsonOutput.toJson(gcta_log.toString())
    pair_log_literal = groovy.json.JsonOutput.toJson(pair_log.toString())
    task_process_literal = groovy.json.JsonOutput.toJson(task.process.toString())
    task_container_literal = groovy.json.JsonOutput.toJson(task.container.toString())
    template('normalise_gcta_bivariate.py')

    stub:
    """
    cat <<-'END_HERITABILITY' > heritability.tsv
    relationship_id\trequest_id\tmethod\tendpoint\tanalysis_id\ttrait_id\ttrait_type\tscale\testimate\tstandard_error\tclassification
    ${meta.relationship_id}\t${meta.request_id}\tgcta_bivariate_reml\tleft\t${meta.left_analysis_id}\t${meta.left_trait_id}\t${meta.left_trait_type}\tobserved\tNA\tNA\tcompleted_nonestimable
    ${meta.relationship_id}\t${meta.request_id}\tgcta_bivariate_reml\tright\t${meta.right_analysis_id}\t${meta.right_trait_id}\t${meta.right_trait_type}\tobserved\tNA\tNA\tcompleted_nonestimable
    END_HERITABILITY
    cat <<-'END_CORRELATION' > genetic_correlation.tsv
    relationship_id\trequest_id\tmethod\tleft_analysis_id\tright_analysis_id\tleft_trait_id\tright_trait_id\testimate\tstandard_error\tclassification
    ${meta.relationship_id}\t${meta.request_id}\tgcta_bivariate_reml\t${meta.left_analysis_id}\t${meta.right_analysis_id}\t${meta.left_trait_id}\t${meta.right_trait_id}\tNA\tNA\tcompleted_nonestimable
    END_CORRELATION
    cat <<-'END_COVARIANCE' > genetic_covariance.tsv
    relationship_id\trequest_id\tmethod\tleft_analysis_id\tright_analysis_id\tleft_trait_id\tright_trait_id\tscale\testimate\tstandard_error\tclassification
    ${meta.relationship_id}\t${meta.request_id}\tgcta_bivariate_reml\t${meta.left_analysis_id}\t${meta.right_analysis_id}\t${meta.left_trait_id}\t${meta.right_trait_id}\tobserved\tNA\tNA\tcompleted_nonestimable
    END_COVARIANCE
    cat <<-'END_DIAGNOSTICS' > diagnostics.tsv
    relationship_id\trequest_id\tmetric\tvalue
    ${meta.relationship_id}\t${meta.request_id}\tclassification\tcompleted_nonestimable
    END_DIAGNOSTICS
    cat <<-'END_PROVENANCE' > provenance.json
    {"schema_version":"1.0","relationship_id":"${meta.relationship_id}","request_id":"${meta.request_id}","method":"gcta_bivariate_reml","classification":"completed_nonestimable","stub":true}
    END_PROVENANCE
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """
}
