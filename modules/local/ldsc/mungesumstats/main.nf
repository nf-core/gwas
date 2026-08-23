process LDSC_MUNGESUMSTATS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(sumstats)
    tuple val(meta2), path(merge_alleles)

    output:
    tuple val(meta), path("${prefix}.sumstats.gz"), emit: munged_sumstats
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

    sumstats_input="${sumstats}"
    if [[ "${sumstats}" == *.gz ]]; then
        gzip --decompress --stdout "${sumstats}" > ldsc_munge_input.tsv
        sumstats_input=ldsc_munge_input.tsv
    fi

    munge_sumstats.py \
        --sumstats "\${sumstats_input}" \
        --merge-alleles "${merge_alleles}" \
        --out "${prefix}" \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    printf 'SNP\tA1\tA2\tN\tZ\n' | gzip -c > "${prefix}.sumstats.gz"
    printf 'LDSC munging stub\n' > "${prefix}.log"
    """
}
