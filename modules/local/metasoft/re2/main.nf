process METASOFT_RE2 {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://ghcr.io/lyh970817/gwas/metasoft@sha256:13c0275d393111daf53ba8e6bfa787627eec6facc7c3b8e4879784497745f073'
        : 'ghcr.io/lyh970817/gwas/metasoft:2.0.1-515de6e@sha256:13c0275d393111daf53ba8e6bfa787627eec6facc7c3b8e4879784497745f073'}"

    input:
    tuple val(meta), path(effect_matrix)

    output:
    tuple val(meta), path("${prefix}.metasoft.txt"), emit: result
    tuple val(meta), path("${prefix}.metasoft.log"), emit: log
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: meta.id
    def heap = (task.memory.toMega() * 0.8).intValue()
    """
    JAVA_TOOL_OPTIONS="-Xmx${heap}M" metasoft \
        ${args} \
        -input "${effect_matrix}" \
        -output "${prefix}.metasoft.txt" \
        -log "${prefix}.metasoft.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metasoft: \$(metasoft --version)
        java: \$(java -version 2>&1 | sed -n 's/.*version "\\([^" ]*\\)".*/\\1/p')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    printf 'RSID\t#STUDY\tPVALUE_FE\tBETA_FE\tSTD_FE\tPVALUE_RE\tBETA_RE\tSTD_RE\tPVALUE_RE2\tSTAT1_RE2\tSTAT2_RE2\tPVALUE_BE\tI_SQUARE\tQ\tPVALUE_Q\tTAU_SQUARE\tPVALUES_OF_STUDIES(Tab_delimitered)\tMVALUES_OF_STUDIES(Tab_delimitered)\\n' > "${prefix}.metasoft.txt"
    printf 'stub\\n' > "${prefix}.metasoft.log"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metasoft: 2.0.1
        java: 17.0.20.1
    END_VERSIONS
    """
}
