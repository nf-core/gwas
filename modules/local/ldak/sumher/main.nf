process LDAK_SUMHER {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://ghcr.io/lyh970817/gwas/ldak@sha256:f2b2157559e4346cab5e9f478ab70fc76359743ef06522fed9ad23769d735a6e'
        : 'ghcr.io/lyh970817/gwas/ldak:6.3-b755ab7@sha256:f2b2157559e4346cab5e9f478ab70fc76359743ef06522fed9ad23769d735a6e'}"

    input:
    tuple val(meta), path(summary_statistics, stageAs: 'summary/*')
    tuple val(meta2), path(tagging_file, stageAs: 'tagging/*')

    output:
    tuple val(meta), path("${prefix}.hers"), emit: hers
    tuple val(meta), path("${prefix}.cats"), emit: categories
    tuple val(meta), path("${prefix}.share"), emit: shares
    tuple val(meta), path("${prefix}.enrich"), emit: enrichments
    tuple val(meta), path("${prefix}.extra"), emit: extra
    tuple val(meta), path("${prefix}.cross"), emit: cross
    tuple val(meta), path("${prefix}.taus"), emit: taus
    tuple val(meta), path("${prefix}.labels"), emit: labels
    tuple val(meta), path("${prefix}.progress"), emit: progress
    tuple val(meta), path("${prefix}.overlap"), emit: overlap
    tuple val(meta), path("${prefix}.hers.liab"), emit: hers_liability, optional: true
    tuple val(meta), path("${prefix}.cats.liab"), emit: categories_liability, optional: true
    tuple val(meta), path("${prefix}.factor"), emit: liability_factor, optional: true
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val("${task.process}"), val('ldak'), eval('cat /opt/ldak/version'), emit: versions_ldak, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: meta.id

    """
    set -o pipefail

    ldak --sum-hers "${prefix}" \
        --summary "${summary_statistics}" \
        --tagfile "${tagging_file}" \
        --max-threads "${task.cpus}" \
        ${args} \
        2>&1 | tee "${prefix}.log"
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: meta.id
    def liability_outputs = args.contains('--prevalence')
        ? """
    touch "${prefix}.hers.liab"
    touch "${prefix}.cats.liab"
    touch "${prefix}.factor"
    """
        : ''

    """
    touch "${prefix}.hers"
    touch "${prefix}.cats"
    touch "${prefix}.share"
    touch "${prefix}.enrich"
    touch "${prefix}.extra"
    touch "${prefix}.cross"
    touch "${prefix}.taus"
    touch "${prefix}.labels"
    touch "${prefix}.progress"
    touch "${prefix}.overlap"
    touch "${prefix}.log"
    ${liability_outputs}
    """
}
