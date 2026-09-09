process GWASLAB_META_ANALYZE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gwaslab:4.1.9--pyhdfd78af_0'
        : 'quay.io/biocontainers/gwaslab:4.1.9--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(sumstats, stageAs: 'parents/*'), val(study_names), val(input_format), val(genome_build), val(random_effects)

    output:
    tuple val(meta), path("${prefix}.fixed.tsv.gz"), emit: fixed
    tuple val(meta), path("${prefix}.random.tsv.gz"), emit: random_effects, optional: true
    tuple val(meta), path("${prefix}.metasoft.input.txt"), emit: metasoft_input
    tuple val(meta), path("${prefix}.mrmega.tsv.gz"), emit: mrmega_input
    tuple val(meta), path("${prefix}.gwaslab.log"), emit: log
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: meta.id
    parents_literal = groovy.json.JsonOutput.toJson(sumstats.collect { parent -> parent.toString() })
    study_names_literal = groovy.json.JsonOutput.toJson(study_names.collect { name -> name.toString() })
    input_format_literal = groovy.json.JsonOutput.toJson(input_format.toString())
    genome_build_literal = groovy.json.JsonOutput.toJson(genome_build.toString())
    random_effects_literal = groovy.json.JsonOutput.toJson(random_effects as boolean)
    prefix_literal = groovy.json.JsonOutput.toJson(prefix.toString())
    template("meta_analyze.py")

    stub:
    prefix = task.ext.prefix ?: meta.id
    def studies = study_names.collect { name -> name.toString() }
    def matrix_fields = (1..studies.size()).collect { _index -> "0 1" }.join(' ')
    def view_header = (1..studies.size()).collect { index -> "EAF_${index}\tBETA_${index}\tSE_${index}\tN_${index}" }.join('\t')
    def random_write = random_effects
        ? "printf 'SNPID\\tCHR\\tPOS\\tEA\\tNEA\\tBETA_RANDOM\\tSE_RANDOM\\tZ_RANDOM\\tP_RANDOM\\n' | gzip -n -c > \"${prefix}.random.tsv.gz\""
        : ''
    """
    printf 'SNPID\tCHR\tPOS\tEA\tNEA\tBETA\tSE\tZ\tP\tDOF\tN\tEAF\tDIRECTION\\n' | gzip -n -c > "${prefix}.fixed.tsv.gz"
    ${random_write}
    printf '1:1000:A:G ${matrix_fields}\\n' > "${prefix}.metasoft.input.txt"
    printf 'META_VARIANT_KEY\tSNPID\tCHR\tPOS\tEA\tNEA\t${view_header}\\n' | gzip -n -c > "${prefix}.mrmega.tsv.gz"
    printf 'stub\\n' > "${prefix}.gwaslab.log"
    cat <<-END_VERSIONS > "versions.yml"
    "${task.process}":
        gwaslab: 4.1.9
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """
}
