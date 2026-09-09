process PREPARE_MRMEGA_INPUT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    tuple val(meta), path(aligned_studies), val(study_names), val(trait_type)

    output:
    tuple val(meta), path("${prefix}.mrmega_inputs/*.txt.gz"), emit: study_files
    tuple val(meta), path("${prefix}.filelist.txt"), emit: filelist
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: meta.id
    aligned_studies_literal = groovy.json.JsonOutput.toJson(aligned_studies.toString())
    study_names_literal = groovy.json.JsonOutput.toJson(study_names.collect { name -> name.toString() })
    trait_type_literal = groovy.json.JsonOutput.toJson(trait_type.toString())
    prefix_literal = groovy.json.JsonOutput.toJson(prefix.toString())
    template("prepare_mrmega_input.py")

    stub:
    prefix = task.ext.prefix ?: meta.id
    def studies = study_names.collect { name -> name.toString() }
    def width = Math.max(2, studies.size().toString().length())
    def effect_header = trait_type.toString() == 'binary'
        ? 'OR\tOR_95L\tOR_95U'
        : 'BETA\tSE'
    def staged = studies
        .withIndex()
        .collect { name, index ->
            "${prefix}.mrmega_inputs/${(index + 1).toString().padLeft(width, '0')}_${name}.txt.gz"
        }
    def study_writes = staged
        .collect { target ->
            "printf 'MARKERNAME\\tEA\\tNEA\\tEAF\\t${effect_header}\\tN\\tCHROMOSOME\\tPOSITION\\n' | gzip -n -c > \"${target}\""
        }
        .join('\n    ')
    def filelist_lines = staged.join('\\n')
    """
    mkdir -p "${prefix}.mrmega_inputs"
    ${study_writes}
    printf '${filelist_lines}\n' > "${prefix}.filelist.txt"
    printf '"${task.process}":\n    python: %s\n' "\$(python3 --version | sed 's/^Python //')" > "versions.yml"
    """
}
