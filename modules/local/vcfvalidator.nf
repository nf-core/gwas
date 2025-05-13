process VCF_VALIDATOR {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::vcftools=0.1.16"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcftools:0.1.16--he513fc3_4' :
        'biocontainers/vcftools:0.1.16--he513fc3_4' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.validator.log"), emit: log
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf-validator \\
        $args \\
        $vcf \\
        > ${prefix}.validator.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.validator.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
    END_VERSIONS
    """
}