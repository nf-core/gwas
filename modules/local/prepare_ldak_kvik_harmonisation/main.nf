// GWASLab requires effect-allele frequency and effective N; the native association table's MAF is not EAF.
// Join A1Freq and n from the native summary table so harmonisation publishes the correct frequency and sample size.
process PREPARE_LDAK_KVIK_HARMONISATION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    // LDAK writes each Step 2 result and its summary table under one shared basename, so the pair is
    // matched by that basename rather than by tuple position.
    tuple val(meta), path(association_files), path(summary_files)

    output:
    tuple val(meta), path('*.harmonisation.tsv'), emit: harmonisation_input
    path 'versions.yml', emit: versions, topic: versions

    script:
    def associations = association_files instanceof List ? association_files : [association_files]
    associations_literal = groovy.json.JsonOutput.toJson(associations.collect { association -> association.toString() })
    task_process_literal = groovy.json.JsonOutput.toJson(task.process.toString())
    template('prepare_ldak_kvik_harmonisation.py')

    stub:
    def associations = association_files instanceof List ? association_files : [association_files]
    def outputs = associations.collect { association -> '"' + association.baseName + '.harmonisation.tsv"' }.join(' ')
    """
    touch ${outputs}
    printf '"${task.process}":\\n    python: %s\\n' "\$(python3 --version | sed 's/^Python //')" > versions.yml
    """
}
