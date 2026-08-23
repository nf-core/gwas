process LDSC_RG {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(left_sumstats, stageAs: 'left/*'), path(right_sumstats, stageAs: 'right/*')
    tuple val(meta2), path(reference_ld_scores, stageAs: 'reference_ld_scores')
    tuple val(meta3), path(regression_weights, stageAs: 'regression_weights')

    output:
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val("${task.process}"), val("ldsc"), eval("python -c 'import importlib.metadata; print(importlib.metadata.version(\"ldsc\"))'"), emit: versions_ldsc, topic: versions
    tuple val("${task.process}"), val("ldsc_native"), eval("cat /opt/venv/ldsc-native-version"), emit: versions_ldsc_native, topic: versions
    tuple val("${task.process}"), val("ldsc_source_revision"), eval("cat /opt/venv/ldsc-source-revision"), emit: versions_ldsc_source, topic: versions
    tuple val("${task.process}"), val("python"), eval("python --version 2>&1 | sed 's/^Python //'"), emit: versions_python, topic: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: meta.id
    """
    export PYTHONUNBUFFERED=1
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export OMP_NUM_THREADS="${task.cpus}"
    export MKL_NUM_THREADS="${task.cpus}"

    ldsc.py \
        --rg "${left_sumstats},${right_sumstats}" \
        --ref-ld-chr "reference_ld_scores/" \
        --w-ld-chr "regression_weights/" \
        --out "${prefix}" \
        ${args} \
        | tee "${prefix}.stdout.log"

    if [[ -s "${prefix}.log" ]]; then
        rm "${prefix}.stdout.log"
    else
        mv "${prefix}.stdout.log" "${prefix}.log"
    fi
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    printf 'LDSC genetic correlation stub\n' > "${prefix}.log"
    """
}
