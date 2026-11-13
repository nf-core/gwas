process CUSTOM_GCTASTRATIFYLDSCORES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/48/483e9d9b3b07e5658792d579e230ad40ed18daf7b9ebfb4323c08570f92fd1d5/data'
        : 'community.wave.seqera.io/library/r-base:4.2.1--b0b5476e2e7a0872'}"

    input:
    tuple val(meta), path(ld_scores), val(ld_bins), val(maf_edges)

    output:
    tuple val(meta), path("${prefix}.strata.tsv"), path("${prefix}_snp_group*.txt"), emit: strata_bundle
    tuple val(meta), path("${prefix}.groups.tsv"), path("${prefix}_snp_group*.txt"), emit: groups
    // `eval()` is the house style for version capture, but Nextflow rejects an `eval` output on any process
    // whose script is not Bash, and this one is an R `template`. The template version and the interpreter
    // version are therefore written from inside the template.
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    task_process_literal = groovy.json.JsonOutput.toJson(task.process.toString())
    template('gcta_stratify_ldscores.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'stub_snp1\\n' > "${prefix}_snp_group_ld01_maf01.txt"
    printf 'stub_snp2\\n' > "${prefix}_snp_group_ld02_maf01.txt"
    printf 'model_key\\tstratum_key\\tld_lower\\tld_upper\\tmaf_lower\\tmaf_upper\\tpredictor_count\\tgroup_filename\\nldms\\tld01_maf01\\t0\\t1\\t0\\t0.5\\t1\\t${prefix}_snp_group_ld01_maf01.txt\\nldms\\tld02_maf01\\t1\\t2\\t0\\t0.5\\t1\\t${prefix}_snp_group_ld02_maf01.txt\\n' > "${prefix}.strata.tsv"
    printf 'group_id\\tgroup_filename\\nld01_maf01\\t${prefix}_snp_group_ld01_maf01.txt\\nld02_maf01\\t${prefix}_snp_group_ld02_maf01.txt\\n' > "${prefix}.groups.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gcta_stratify_ldscores: 1.0.0
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """
}
