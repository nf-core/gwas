process CUSTOM_MPHSNPINFO {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
        : 'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927'}"

    input:
    tuple val(meta), path(bim)
    // The plan is one metadata-bearing tuple because its manifest and its group files are only meaningful
    // together. `weight_names` is the ordered column list the caller wants and is never derived here: it is a
    // mandatory tuple-local scalar, echoed back on the output tuple so the consumer reads the column order
    // from the contract rather than by re-parsing the CSV. `[]`, `[]` with `['all']` is the one-component case.
    tuple val(meta2), path(strata_manifest), path(snp_group_files), val(weight_names)
    // Process-global, the shape `GCTA_CALCULATELDSCORES` already uses for its region size. It is the highest
    // integer chromosome code that counts as autosomal, matching GCTA's `--autosome-num`.
    val autosome_count

    output:
    // Two emits rather than one bundle: the CSV is what MPH reads through `--snp_info_file`, while the counts
    // are pipeline accounting that no native command ever sees, and no consumer wants them together.
    tuple val(meta), path("${prefix}.snp_info.csv"), val(weight_names), emit: snp_info
    tuple val(meta), path("${prefix}.snp_info.counts.tsv"), emit: counts
    tuple val("${task.process}"), val("python"), eval("python3 --version | sed 's/^Python //'"), emit: versions_python, topic: versions


    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def manifest_argument = strata_manifest ? strata_manifest.toString() : ''
    """
    python3 - \\
        "${bim}" \\
        "${prefix}" \\
        "${autosome_count}" \\
        "${manifest_argument}" \\
        "${weight_names.join(',')}" <<'PY'
    import csv
    import sys

    bim_path, prefix, autosome_count, manifest_path, weight_name_list = sys.argv[1:6]
    autosome_count = int(autosome_count)
    columns = [name for name in weight_name_list.split(",") if name != ""]


    def fail(message):
        sys.exit("[nf-core/gwas] ERROR: {}".format(message))


    if not columns:
        fail("CUSTOM_MPHSNPINFO was given no weight column names")
    if len(set(columns)) != len(columns):
        fail("weight column names repeat: {}".format(columns))

    variants = []
    with open(bim_path) as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            variants.append((fields[0], fields[1]))

    # The universe is declared rather than inherited from the BIM. GCTA restricts `--make-grm` to the autosomes
    # and refuses `--ld-score-region` outright on anything else, while MPH applies no chromosome rule at all and
    # says nothing about it: on one bfile carrying 100 non-autosomal variants GCTA used 2100 and MPH used 2200,
    # both at exit 0. A variant outside the universe is written with weight 0 in every column, which is the same
    # representation as "in no stratum", so the CSV stays auditable against the BIM line for line.
    def in_universe(chromosome):
        try:
            code = int(chromosome)
        except ValueError:
            return False
        return 1 <= code <= autosome_count

    known_variants = set(name for _chromosome, name in variants)
    membership = {}

    if manifest_path:
        with open(manifest_path) as handle:
            manifest_rows = [row for row in handle.read().splitlines() if row != ""]
        header = manifest_rows[0].split("\\t")
        records = [dict(zip(header, row.split("\\t"))) for row in manifest_rows[1:]]
        declared = [record["stratum_key"] for record in records]
        if declared != columns:
            fail(
                "the caller's weight column order {} does not match the plan manifest's stratum order {}".format(
                    columns, declared
                )
            )
        universe_by_name = {name: in_universe(chromosome) for chromosome, name in variants}
        for index, record in enumerate(records):
            with open(record["group_filename"]) as handle:
                group = [name for name in handle.read().split() if name != ""]
            for name in group:
                if name not in known_variants:
                    fail(
                        "SNP '{}' is in group file '{}' but not in '{}'; the component plan and the genotype "
                        "bundle are not the same view".format(name, record["group_filename"], bim_path)
                    )
                if not universe_by_name[name]:
                    fail(
                        "SNP '{}' is in group file '{}' but lies outside the autosomal universe (chromosome "
                        "codes 1 to {})".format(name, record["group_filename"], autosome_count)
                    )
                if name in membership:
                    fail(
                        "SNP '{}' is in more than one group file ('{}' and '{}'); the component plan declares "
                        "disjoint groups".format(name, records[membership[name]]["group_filename"], record["group_filename"])
                    )
                membership[name] = index
    else:
        if columns != ["all"]:
            fail("without a plan manifest exactly one weight column named 'all' is expected, got {}".format(columns))
        for chromosome, name in variants:
            if in_universe(chromosome):
                membership[name] = 0

    counts = [0] * len(columns)
    unassigned = 0
    with open("{}.snp_info.csv".format(prefix), "w", newline="\\n") as handle:
        writer = csv.writer(handle, lineterminator="\\n")
        writer.writerow(["SNP"] + columns)
        for _chromosome, name in variants:
            index = membership.get(name)
            if index is None:
                unassigned += 1
            else:
                counts[index] += 1
            writer.writerow([name] + ["1" if index == column else "0" for column in range(len(columns))])

    # `__unassigned__` is the number of BIM rows MPH is told to ignore and `__bim_rows__` the number it will
    # report as found, so the pair reconciles the pipeline's declared predictor set against MPH's own echo.
    with open("{}.snp_info.counts.tsv".format(prefix), "w", newline="\\n") as handle:
        handle.write("stratum_key\\tpredictor_count\\n")
        for name, count in zip(columns, counts):
            handle.write("{}\\t{}\\n".format(name, count))
        handle.write("__unassigned__\\t{}\\n".format(unassigned))
        handle.write("__bim_rows__\\t{}\\n".format(len(variants)))
    PY
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // The stub's column count follows the caller's `weight_names`, so a stub run of a K-component plan yields
    // a K-column CSV and the calling composition's stub topology is exercised for real.
    def stub_header = (['SNP'] + weight_names).join(',')
    def stub_row = (['stub_snp1'] + weight_names.withIndex().collect { _name, index -> index == 0 ? '1' : '0' }).join(',')
    def stub_counts = weight_names.withIndex().collect { name, index -> "${name}\\t${index == 0 ? 1 : 0}" }.join('\\n')
    """
    printf '%s\\n%s\\n' '${stub_header}' '${stub_row}' > "${prefix}.snp_info.csv"
    printf 'stratum_key\\tpredictor_count\\n${stub_counts}\\n__unassigned__\\t0\\n__bim_rows__\\t1\\n' > "${prefix}.snp_info.counts.tsv"
    """
}
