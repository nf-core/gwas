process CUSTOM_MPHSNPINFO {
    tag "${meta.id}_${meta2.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    tuple val(meta), path(bim)
    // Named group files and their manifest are staged together; row order is the weight-column order.
    tuple val(meta2), path(group_manifest), path(snp_group_files), val(weight_names)

    output:
    // Two emits rather than one bundle: the CSV is what MPH reads through `--snp_info_file`, while the counts
    // are group accounting that no native command ever sees, and no consumer wants them together.
    tuple val(meta), path("${prefix}.snp_info.csv"), val(weight_names), emit: snp_info
    tuple val(meta), path("${prefix}.snp_info.counts.tsv"), emit: counts
    // `eval()` is the house style for version capture, but Nextflow rejects an `eval` output on any process
    // whose script is not Bash, and this one is a Python `template`. The interpreter version is therefore
    // written from inside the template.
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    bim_literal = groovy.json.JsonOutput.toJson(bim.toString())
    prefix_literal = groovy.json.JsonOutput.toJson(prefix.toString())
    manifest_literal = groovy.json.JsonOutput.toJson(group_manifest.toString())
    columns_literal = groovy.json.JsonOutput.toJson(weight_names)
    task_process_literal = groovy.json.JsonOutput.toJson(task.process.toString())
    template('mphsnpinfo.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // The stub's column count follows the caller's `weight_names`, so a stub run of a K-component plan yields
    // a K-column CSV and the calling composition's stub topology is exercised for real.
    def stub_header = (['SNP'] + weight_names).join(',')
    def stub_row = (['stub_snp1'] + weight_names.withIndex().collect { _name, index -> index == 0 ? '1' : '0' }).join(',')
    def stub_counts = weight_names.withIndex().collect { name, index -> "${name}\\t${index == 0 ? 1 : 0}" }.join('\\n')
    """
    printf '%s\\n%s\\n' '${stub_header}' '${stub_row}' > "${prefix}.snp_info.csv"
    printf 'group_id\\tpredictor_count\\n${stub_counts}\\n__unassigned__\\t0\\n__bim_rows__\\t1\\n' > "${prefix}.snp_info.counts.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """
}
