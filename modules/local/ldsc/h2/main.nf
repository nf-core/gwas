process LDSC_H2 {
    tag "${meta.id}"
    label 'process_low'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://ghcr.io/lyh970817/gwas/ldsc@sha256:77fbb697c16a559c3fe75204b1e7ab6a0202afcf10b8a8629bcc98592b0e412b'
        : 'ghcr.io/lyh970817/gwas/ldsc:3.0.2-cbiit-6c67395@sha256:77fbb697c16a559c3fe75204b1e7ab6a0202afcf10b8a8629bcc98592b0e412b'}"

    input:
    tuple val(meta), path(sumstats)
    tuple val(meta2), path(reference_ld_scores, stageAs: 'reference_ld_scores')
    tuple val(meta3), path(regression_weights, stageAs: 'regression_weights')

    output:
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val("${task.process}"), val("ldsc"), eval("cat /opt/venv/ldsc-source-revision"), emit: versions_ldsc, topic: versions
    tuple val("${task.process}"), val("python"), eval("python --version 2>&1 | sed 's/^Python //'"), emit: versions_python, topic: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: meta.id
    """
    export PYTHONUNBUFFERED=1
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export OMP_NUM_THREADS="${task.cpus}"
    export MKL_NUM_THREADS="${task.cpus}"

    # This LDSC release writes the complete native analysis log to stdout. Retain that stream under the
    # declared result name; the file it opens at the same name is empty. Retire when the pinned CBIIT ldsc39
    # revision moves to one whose --out .log carries the estimates.
    ldsc.py \
        --h2 "${sumstats}" \
        --ref-ld-chr "reference_ld_scores/" \
        --w-ld-chr "regression_weights/" \
        --out "${prefix}" \
        ${args} \
        > "${prefix}.stdout.log"

    mv "${prefix}.stdout.log" "${prefix}.log"
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    printf 'LDSC heritability stub\n' > "${prefix}.log"
    """
}
