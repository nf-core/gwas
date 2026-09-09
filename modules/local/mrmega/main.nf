process MRMEGA {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://quay.io/loukas_moutsianas/mrmega@sha256:1143b7f016f00f0f32cbc6ad72b4a571e2f824e440a0b7fecf4ec4c8334ae1e8'
        : 'quay.io/loukas_moutsianas/mrmega@sha256:1143b7f016f00f0f32cbc6ad72b4a571e2f824e440a0b7fecf4ec4c8334ae1e8'}"

    input:
    tuple val(meta), path(study_files, stageAs: 'studies/*'), path(filelist), val(axes), val(trait_type)

    output:
    tuple val(meta), path("${prefix}.result"), emit: result
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: meta.id
    def quantitative_flag = trait_type.toString() == 'quantitative' ? '--qt' : ''
    """
    : > "${prefix}.mrmega.in"
    while IFS= read -r source; do
        printf 'studies/%s\n' "\$(basename "\$source")" >> "${prefix}.mrmega.in"
    done < "${filelist}"

    /MR-MEGA/MR-MEGA \
        -i "${prefix}.mrmega.in" \
        --pc "${axes}" \
        ${quantitative_flag} \
        -o "${prefix}"

    printf '"%s":\n    mrmega: %s\n' \
        "${task.process}" \
        "\$(/MR-MEGA/MR-MEGA --version 2>&1 | sed -n 's/.*version: *//p' | head -n 1)" \
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    def stub_axes = axes as int
    def coefficients = (0..stub_axes).collect { index -> "beta_${index}\tse_${index}" }.join('\t')
    def stub_header = [
        'MarkerName\tChromosome\tPosition\tEA\tNEA\tEAF\tNsample\tNcohort\tEffects',
        coefficients,
        'chisq_association\tndf_association\tP-value_association',
        'chisq_ancestry_het\tndf_ancestry_het\tP-value_ancestry_het',
        'chisq_residual_het\tndf_residual_het\tP-value_residual_het',
        'lnBF\tComments',
    ].join('\t')
    """
    printf '%s\n' '${stub_header}' > "${prefix}.result"
    printf '###################\n# MR-MEGA v.0.2\n###################\n\nAnalysis finished.\n' > "${prefix}.log"
    printf '"%s":\n    mrmega: 0.2\n' "${task.process}" > versions.yml
    """
}
