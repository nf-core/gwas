process FINALISE_META_ANALYSIS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gwaslab:4.1.9--pyhdfd78af_0'
        : 'quay.io/biocontainers/gwaslab:4.1.9--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fixed), path(random_result), path(metasoft), path(mrmega)

    output:
    tuple val(meta), path("${prefix}.gwaslab.tsv.gz"), emit: summary_statistics
    path 'versions.yml', emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: meta.summary_statistics_id
    fixed_literal = groovy.json.JsonOutput.toJson(fixed.toString())
    random_literal = groovy.json.JsonOutput.toJson(random_result ? random_result.toString() : null)
    metasoft_literal = groovy.json.JsonOutput.toJson(metasoft ? metasoft.toString() : null)
    mrmega_literal = groovy.json.JsonOutput.toJson(mrmega ? mrmega.toString() : null)
    min_studies = meta.min_studies
    output_literal = groovy.json.JsonOutput.toJson("${prefix}.gwaslab.tsv.gz")
    task_process_literal = groovy.json.JsonOutput.toJson(task.process.toString())
    template('finalise_meta_analysis.py')

    stub:
    prefix = task.ext.prefix ?: meta.summary_statistics_id
    def optional_columns = []
    if (random_result) {
        optional_columns.addAll(['BETA_RANDOM', 'SE_RANDOM', 'Z_RANDOM', 'P_RANDOM'])
    }
    if (metasoft) {
        optional_columns.addAll(['P_RE2', 'RE2_MEAN_COMPONENT', 'RE2_HET_COMPONENT'])
    }
    if (mrmega) {
        optional_columns.addAll(['MRMEGA_CHISQ_ASSOC', 'MRMEGA_DF_ASSOC', 'MRMEGA_P_ASSOC', 'MRMEGA_CHISQ_ANCESTRY_HET', 'MRMEGA_DF_ANCESTRY_HET', 'MRMEGA_P_ANCESTRY_HET', 'MRMEGA_CHISQ_RESIDUAL_HET', 'MRMEGA_DF_RESIDUAL_HET', 'MRMEGA_P_RESIDUAL_HET', 'MRMEGA_LNBF'])
    }
    def header = (['SNPID', 'CHR', 'POS', 'EA', 'NEA', 'BETA', 'SE', 'P', 'EAF', 'N', 'Z', 'Q', 'P_HET', 'I2', 'DIRECTION', 'DOF', 'N_STUDIES'] + optional_columns).join('\\t')
    """
    printf '${header}\\n' | gzip -n -c > "${prefix}.gwaslab.tsv.gz"
    python3 -c 'import platform, pandas; print("\\\"${task.process}\\\":\\n    python: " + platform.python_version() + "\\n    pandas: " + pandas.__version__)' > versions.yml
    """
}
