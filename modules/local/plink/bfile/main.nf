    process PLINK_BFILE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::plink=1.90b6.21" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1':
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(bed), path(bim),path(fam)
    path(dups)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.bim"), emit: bim
    tuple val(meta), path("*.fam"), emit: fam
    tuple val(meta), path("*.nosex"), emit: nosex, optional: true
    tuple val(meta), path("*.dups"), emit: dups, optional: true
    tuple val(meta), path("*.orig"), emit: orig, optional: true
    tuple val(meta), path("*.imiss"), emit: imiss, optional: true
    tuple val(meta), path("*.hwe"), emit: hwe, optional: true
    tuple val(meta), path("*.het"), emit: het, optional: true
    tuple val(meta), path("*.freq"), emit: freq, optional: true
    tuple val(meta), path("*.irem"), emit: irem, optional: true
    tuple val(meta), path("*.eigenval"), emit: eigenval, optional: true
    tuple val(meta), path("*.eigenvec"), emit: eigenvac, optional: true
    tuple val(meta), path("*.genome"), emit: genome, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def exlude_dups = "${dups}" ? "--exclude ${dups}" : ""
    """
    plink \\
        --bfile ${meta.id} \\
        --threads $task.cpus \\
        ${exclude_dups} \\
        ${args} \\
        --out ${prefix}_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """
}
