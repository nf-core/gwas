process MPH_MAKEGRM {
    tag "${meta.id}_${meta2.id}"
    label 'process_high'

    // MPH is not packaged for Bioconda, so this module ships a digest-pinned image and no `environment.yml`,
    // following `metasoft/re2`. The image is `linux/amd64` only, statically linked against Intel
    // oneMKL, and there is no Conda fallback, so the module cannot run on another architecture.
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://ghcr.io/lyh970817/gwas/mph@sha256:c471ae7132e961bf401dd3db177d4c0c1c6cfdf0fc8306db1dcaac5b85a6c3e7'
        : 'ghcr.io/lyh970817/gwas/mph:0.55.1-13ffe63@sha256:c471ae7132e961bf401dd3db177d4c0c1c6cfdf0fc8306db1dcaac5b85a6c3e7'}"

    input:
    // `bed`, `bim` and `fam` share one basename: MPH takes a `--bfile` prefix, not three paths. That is a
    // caller contract stated in `meta.yml` rather than a defensive check here.
    tuple val(meta), path(bed), path(bim), path(fam)
    // `--snp_weight_name` is mandatory natively: MPH exits 1 naming the column when it is absent from the
    // header, and builds nothing when the flag is omitted. It is therefore a tuple-local scalar, not `ext.args`.
    tuple val(meta2), path(snp_info), val(weight_name)

    output:
    // `.grm.bin` and `.grm.iid` are the whole MPH matrix and are never consumed apart, so they travel as one
    // collected path. There is no `.grm.N.bin`: the binary's header carries the sum of the SNP weights, and
    // the triangle below it is the raw unnormalised cross-product.
    tuple val(meta), path("${prefix}.grm.{bin,iid}"), emit: grm_files
    tuple val(meta), path("${prefix}.log"), emit: log
    // MPH has no --version; the bare binary prints its banner and exits 1, tolerated by || true.
    tuple val("${task.process}"), val("mph"), eval("(mph 2>&1 || true) | sed -n 's/^[*] Version \\([0-9][0-9.]*\\).*/\\1/p'"), emit: versions_mph, topic: versions

    script:
    def args = task.ext.args ?: ''
    def bfile_prefix = bed.baseName
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mph \\
        --make_grm \\
        --bfile "${bfile_prefix}" \\
        --snp_info_file "${snp_info}" \\
        --snp_weight_name "${weight_name}" \\
        --output_file "${prefix}" \\
        --num_threads "${task.cpus}" \\
        ${args} \\
        2>&1 | tee "${prefix}.log"

    # Retire when an MPH release returns a non-zero
    # exit on failure.
    # MPH's `main()` catches most thrown errors, prints them and returns 0, so the exit status
    # is not the contract: the primary bundle and the absence of an error line in the log are. Both members are asserted
    # because a missing `.grm.iid` alone would make every downstream sample map silently empty.
    test -s "${prefix}.grm.bin"
    test -s "${prefix}.grm.iid"
    ! grep -qE '^(Error|Inconsistency)' "${prefix}.log"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // A parseable header rather than an empty file: `int32 n = 200` followed by a `float32` weight sum of 0,
    // which is the layout every consumer of this bundle reads first.
    """
    printf '\\310\\000\\000\\000\\000\\000\\000\\000' > "${prefix}.grm.bin"
    printf 'stub_iid\\n' > "${prefix}.grm.iid"
    printf 'stub\\n' > "${prefix}.log"
    """
}
