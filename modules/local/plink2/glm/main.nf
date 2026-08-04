process PLINK2_GLM {
    tag "${meta.id}_${meta2.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/22/22dd30b4bc53d1ba88fcd146a2dd562659c9aaa3f0f1712542fb2b4048b454b7/data'
        : 'community.wave.seqera.io/library/plink2:2.0.0a.6.9--e6710830a4b7f0c6'}"

    input:
    tuple val(meta), path(pgen), path(psam), path(pvar)
    tuple val(meta2), path(phenotype)
    tuple val(meta3), path(covariates)

    output:
    tuple val(meta), path("${prefix}.*.glm.linear"), emit: linear, optional: true
    tuple val(meta), path("${prefix}.*.glm.logistic"), emit: logistic, optional: true
    tuple val(meta), path("${prefix}.*.glm.logistic.hybrid"), emit: logistic_hybrid, optional: true
    tuple val(meta), path("${prefix}.*.glm.firth"), emit: firth, optional: true
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val("${task.process}"), val("plink2"), eval("plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//'"), emit: versions_plink2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_prefix = pgen.baseName
    prefix = task.ext.prefix ?: "${meta.id}"
    def covariates_arg = covariates ? "--covar \"${covariates}\"" : ''
    def mem_mb = task.memory.toMega()
    """
    plink2 \\
        --pfile "${input_prefix}" \\
        --pheno "${phenotype}" \\
        ${covariates_arg} \\
        --threads "${task.cpus}" \\
        --memory "${mem_mb}" \\
        --out "${prefix}" \\
        --glm \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP" > "${prefix}.stub.glm.linear"
    echo "#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tOR\tLOG(OR)_SE\tZ_STAT\tP" > "${prefix}.stub.glm.logistic"
    echo "#CHROM\tPOS\tID\tREF\tALT\tA1\tFIRTH?\tTEST\tOBS_CT\tOR\tLOG(OR)_SE\tZ_STAT\tP" > "${prefix}.stub.glm.logistic.hybrid"
    echo "#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tOR\tLOG(OR)_SE\tZ_STAT\tP" > "${prefix}.stub.glm.firth"
    echo "PLINK 2 GLM stub" > "${prefix}.log"
    """
}
