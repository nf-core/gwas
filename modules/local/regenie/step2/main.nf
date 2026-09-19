process REGENIE_STEP2 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7a/7a05bf71ea09adc5ebf9f0c656c9b326c0f16ba8e4966914972e58313469a466/data'
        : 'community.wave.seqera.io/library/regenie:4.1.2--5d361f9fcb2f85cf'}"

    input:
    tuple val(meta), path(plink_genotype_file), path(plink_variant_file), path(plink_sample_file)
    tuple val(meta2), path(predictions), path(loco_files)
    tuple val(meta3), path(pheno)
    tuple val(meta4), path(covar)
    val bsize

    output:
    tuple val(meta), path("*.regenie.gz"), emit: results
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('regenie'), eval('regenie --version 2>&1 | sed -n "1{s/^v//;s/\\.gz$//;p}"'), topic: versions, emit: versions_regenie

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def covar_arg = covar ? "--covarFile ${covar}" : ''
    def genotype_flag = plink_genotype_file.name.endsWith('.pgen') ? '--pgen' : '--bed'
    def input_prefix = plink_genotype_file.baseName
    def prefix = task.ext.prefix ?: "${meta.id}.${plink_genotype_file.baseName}"

    """
    regenie \\
        --step 2 \\
        ${genotype_flag} ${input_prefix} \\
        --phenoFile ${pheno} \\
        --pred ${predictions} \\
        ${covar_arg} \\
        --bsize ${bsize} \\
        --gz \\
        --threads ${task.cpus} \\
        ${args} \\
        --out ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${plink_genotype_file.baseName}"
    """
    index=0
    for phenotype in \$(awk 'NR == 1 { for (i = 3; i <= NF; i++) print \$i; exit }' ${pheno}); do
        index=\$((index + 1))
        echo "" | gzip > ${prefix}_\${phenotype}.regenie.gz
    done
    touch ${prefix}.log
    """
}
