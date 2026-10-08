process GWASLAB_HARMONIZE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gwaslab:4.1.9--pyhdfd78af_0'
        : 'quay.io/biocontainers/gwaslab:4.1.9--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(sumstats), val(input_format), val(genome_build)
    tuple val(meta2), path(reference_fasta), path(reference_fasta_fai)
    tuple val(meta3), path(rsid_reference_vcf), path(rsid_reference_vcf_index)
    tuple val(meta4), path(strand_reference_vcf), path(strand_reference_vcf_index)

    output:
    tuple val(meta), path("*.gwaslab.tsv.gz"), emit: sumstats
    tuple val(meta), path("*.gwaslab.log"), emit: log
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template("harmonize.py")

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '' | gzip -c > "${prefix}.gwaslab.tsv.gz"
    touch "${prefix}.gwaslab.log"
    cat <<-END_VERSIONS > "versions.yml"
    "${task.process}":
        gwaslab: 4.1.9
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """
}
