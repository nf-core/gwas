process LDAK_SUMCORS {
    tag "${meta.id}_${meta2.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://ghcr.io/lyh970817/gwas/ldak@sha256:a0699bb7637ef2ba20c37fc6ef46b9533c78e68a022d770478a2fbb938e10a00'
        : 'ghcr.io/lyh970817/gwas/ldak:6.3-b755ab7-bookworm@sha256:a0699bb7637ef2ba20c37fc6ef46b9533c78e68a022d770478a2fbb938e10a00'}"

    input:
    tuple val(meta), path(summary_statistics, stageAs: 'trait1/*')
    tuple val(meta2), path(summary_statistics2, stageAs: 'trait2/*')
    tuple val(meta3), path(tagging_file, stageAs: 'tagging/*')

    output:
    tuple val(meta), val(meta2), path("${prefix}.cors"), emit: correlations
    tuple val(meta), val(meta2), path("${prefix}.cors.full"), emit: correlations_full
    tuple val(meta), val(meta2), path("${prefix}.labels"), emit: labels
    tuple val(meta), val(meta2), path("${prefix}.progress"), emit: progress
    tuple val(meta), val(meta2), path("${prefix}.overlap"), emit: overlap
    tuple val(meta), val(meta2), path("${prefix}.cors.liab"), emit: correlations_liability, optional: true
    tuple val(meta), val(meta2), path("${prefix}.log"), emit: log
    tuple val("${task.process}"), val('ldak'), eval('cat /opt/ldak/version'), emit: versions_ldak, topic: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_${meta2.id}"

    """
    ldak --sum-cors "${prefix}" \
        --summary "${summary_statistics}" \
        --summary2 "${summary_statistics2}" \
        --tagfile "${tagging_file}" \
        --max-threads "${task.cpus}" \
        ${args} \
        2>&1 | tee "${prefix}.log"
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_${meta2.id}"
    def liability_outputs = args.contains('--prevalence')
        ? """
    touch "${prefix}.cors.liab"
    """
        : ''

    """
    touch "${prefix}.cors"
    touch "${prefix}.cors.full"
    touch "${prefix}.labels"
    touch "${prefix}.progress"
    touch "${prefix}.overlap"
    touch "${prefix}.log"
    ${liability_outputs}
    """
}
