process GCTA_MAKEGRMPART {
    tag "${meta.id}: part ${part} of ${nparts}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46b0d05f0daa47561d87d2a9cac5e51edc2c78e26f1bbab439c688386241a274/data'
        : 'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9'}"

    input:
    tuple val(meta), path(mfile), path(bed_pgen), path(bim_pvar), path(fam_psam), val(nparts), val(part)
    tuple val(meta2), path(snp_group_file)

    output:
    // GCTA zero-pads the part index to the digit width of the requested part count, so `--make-grm-part 100 3`
    // writes `.part_100_003.grm.*` while `--make-grm-part 9 1` writes `.part_9_1.grm.*`. The index is emitted
    // as a value beside the files, so the glob does not have to restate the padding rule.
    tuple val(meta), path("*.part_${nparts}_*.grm.*"), val(nparts), val(part), emit: grm_files
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val("gcta"), eval("gcta --version | sed -En 's/^[*] version v([0-9.]*).*/\\1/p'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genotype_files = bed_pgen instanceof List ? bed_pgen : [bed_pgen]
    def genotype_extension = genotype_files[0].name.tokenize('.').last()
    def prefix = task.ext.prefix ?: "${meta.id}"
    def multi_file_flag = genotype_extension == 'pgen' ? '--mpfile' : '--mbfile'
    def extract_cmd = snp_group_file ? "--extract \"${snp_group_file}\"" : ''
    """
    gcta \\
        ${multi_file_flag} "${mfile}" \\
        --make-grm-part "${nparts}" "${part}" \\
        ${extract_cmd} \\
        --thread-num "${task.cpus}" \\
        --out "${prefix}" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def part_index = String.format("%0${nparts.toString().length()}d", part as int)
    """
    touch "${prefix}.part_${nparts}_${part_index}.grm.id"
    touch "${prefix}.part_${nparts}_${part_index}.grm.bin"
    touch "${prefix}.part_${nparts}_${part_index}.grm.N.bin"
    touch "${prefix}.part_${nparts}_${part_index}.log"
    """
}
