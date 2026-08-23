process CANONICALISE_SUMMARY_STATISTICS {
    tag "${meta.summary_statistics_id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    tuple val(meta), path(candidate, stageAs: 'candidate/*'), path(source, stageAs: 'source/*')

    output:
    tuple val(meta), path("${prefix}.canonical.tsv.gz"), emit: summary_statistics
    tuple val(meta), path("${prefix}.provenance.json"), emit: provenance
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: meta.summary_statistics_id
    candidate_literal = groovy.json.JsonOutput.toJson(candidate.toString())
    source_literal = groovy.json.JsonOutput.toJson(source.toString())
    metadata_literal = groovy.json.JsonOutput.toJson(groovy.json.JsonOutput.toJson(meta))
    output_literal = groovy.json.JsonOutput.toJson("${prefix}.canonical.tsv.gz")
    provenance_literal = groovy.json.JsonOutput.toJson("${prefix}.provenance.json")
    template("canonicalise.py")

    stub:
    prefix = task.ext.prefix ?: meta.summary_statistics_id
    def metadata = groovy.json.JsonOutput.toJson(
        [
            summary_statistics_id: meta.summary_statistics_id,
            canonical_contract_version: 'nfcore_gwas_canonical_v1',
            stub: true,
        ]
    )
    """
    printf 'SNPID\tCHR\tPOS\tEA\tNEA\tSTATUS\tEAF\tBETA\tSE\tP\tN\\n' | gzip -n -c > "${prefix}.canonical.tsv.gz"
    printf '%s\\n' '${metadata}' > "${prefix}.provenance.json"
    printf '"${task.process}":\\n    python: %s\\n' "\$(python3 --version | sed 's/^Python //')" > "versions.yml"
    """
}
