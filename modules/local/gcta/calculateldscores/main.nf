process GCTA_CALCULATELDSCORES {
    tag "${chr ? "${meta.id} chr${chr}" : meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46b0d05f0daa47561d87d2a9cac5e51edc2c78e26f1bbab439c688386241a274/data'
        : 'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9'}"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), val(ld_score_region), val(chr)

    output:
    tuple val(meta), path("*_gcta_ld.score.ld"), emit: ld_scores
    tuple val(meta), path("*_gcta_ld.log"), emit: log
    tuple val("${task.process}"), val("gcta"), eval("gcta --version | sed -En 's/^[*] version v([0-9.]*).*/\\1/p'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def bfile_prefix = bed.baseName
    // A per-chromosome call gets its own default prefix, so the tables of one bundle's chromosomes can be
    // staged side by side.
    def prefix = task.ext.prefix ?: (chr ? "${meta.id}_chr${chr}" : "${meta.id}")
    def chr_arg = chr ? "--chr ${chr}" : ''
    """
    gcta \\
        --bfile "${bfile_prefix}" \\
        ${chr_arg} \\
        --ld-score-region "${ld_score_region}" \\
        --out "${prefix}_gcta_ld" \\
        --thread-num "${task.cpus}" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: (chr ? "${meta.id}_chr${chr}" : "${meta.id}")
    """
    printf "%s\\n" \
        "SNP chr bp freq mean_rsq snp_num max_rsq ldscore_SNP ldscore_region" \
        "stub_snp1 1 100 0.10 0.01 10 0.20 1.10 1.20" \
        "stub_snp2 1 200 0.20 0.02 10 0.30 1.20 1.30" \
        > "${prefix}_gcta_ld.score.ld"
    touch "${prefix}_gcta_ld.log"
    """
}
