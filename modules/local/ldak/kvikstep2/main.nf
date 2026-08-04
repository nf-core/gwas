process LDAK_KVIKSTEP2 {
    tag "${meta.id}:${meta2.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c94e5424f08cbd7e0856ab3c3a9992b4080944a5d0c497ce0abcceb413db4a3e/data'
        : 'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129'}"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(phenotype_file)
    tuple val(meta3), path(step1_root), path(step1_loco_details), path(step1_loco_prs)
    tuple val(meta4), path(quant_covariates_file)
    tuple val(meta5), path(cat_covariates_file)
    tuple val(meta6), path(keep_file)

    output:
    tuple val(meta), path("${prefix}.step2*.assoc"), emit: results
    tuple val(meta), path("${prefix}.step2*.harmonisation.tsv"), emit: harmonisation_input
    tuple val(meta), path("${prefix}.step2*.summaries"), emit: summaries
    tuple val(meta), path("${prefix}.step2*.pvalues"), emit: pvalues, optional: true
    tuple val(meta), path("${prefix}.step2.log"), emit: log
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | grep -oP '(?<=^Version )[0-9.]+'"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def bfile_prefix = bed.baseName
    prefix = task.ext.prefix ?: meta2.id
    def covar_arg = quant_covariates_file ? "--covar \"${quant_covariates_file}\"" : ""
    def factors_arg = cat_covariates_file ? "--factors \"${cat_covariates_file}\"" : ""
    def keep_arg = keep_file ? "--keep \"${keep_file}\"" : ""
    """
    if [[ "${step1_root}" != "${prefix}.step1.root" ]]; then
        ln -s "${step1_root}" "${prefix}.step1.root"
        ln -s "${step1_loco_details}" "${prefix}.step1.loco.details"
        ln -s "${step1_loco_prs}" "${prefix}.step1.loco.prs"
    fi

    ldak6 --kvik-step2 "${prefix}" \\
        --bfile "${bfile_prefix}" \\
        --pheno "${phenotype_file}" \\
        ${covar_arg} \\
        ${factors_arg} \\
        ${keep_arg} \\
        --max-threads "${task.cpus}" \\
        ${args} \\
        2>&1 | tee "${prefix}.step2.log"

    for assoc in "${prefix}".step2*.assoc; do
        summary="\${assoc%.assoc}.summaries"
        harmonisation_input="\${assoc%.assoc}.harmonisation.tsv"
        if [[ ! -f "\${summary}" ]]; then
            echo "LDAK-KVIK did not write the summary table required to recover per-variant effect-allele frequencies and effective analysis sizes: \${summary}" >&2
            exit 1
        fi
        awk '
            BEGIN { FS = OFS = "\\t" }
            NR == FNR {
                if (FNR == 1) {
                    for (i = 1; i <= NF; i++) {
                        if (\$i == "Predictor") summary_predictor = i
                        if (\$i == "A1") summary_a1 = i
                        if (\$i == "A2") summary_a2 = i
                        if (\$i == "n") summary_n = i
                        if (\$i == "A1Freq") summary_a1freq = i
                    }
                    if (!summary_predictor || !summary_a1 || !summary_a2 || !summary_n || !summary_a1freq) {
                        print "LDAK-KVIK summary table lacks Predictor, A1, A2, n or A1Freq" > "/dev/stderr"
                        exit 1
                    }
                    next
                }
                key = \$summary_predictor SUBSEP \$summary_a1 SUBSEP \$summary_a2
                if (key in sample_count) {
                    print "Duplicate Predictor+A1+A2 key in LDAK-KVIK summary table: " \$summary_predictor, \$summary_a1, \$summary_a2 > "/dev/stderr"
                    exit 1
                }
                sample_count[key] = \$summary_n
                effect_allele_frequency[key] = \$summary_a1freq
                next
            }
            FNR == 1 {
                for (i = 1; i <= NF; i++) {
                    header[i] = \$i
                    if (\$i == "Predictor") association_predictor = i
                    if (\$i == "A1") association_a1 = i
                    if (\$i == "A2") association_a2 = i
                }
                if (!association_predictor || !association_a1 || !association_a2) {
                    print "LDAK-KVIK association table lacks Predictor, A1 or A2" > "/dev/stderr"
                    exit 1
                }
                first = 1
                for (i = 1; i <= NF; i++) {
                    if (\$i != "Wald_Stat" && \$i != "MAF") {
                        printf "%s%s", (first ? "" : OFS), \$i
                        first = 0
                    }
                }
                print OFS "EAF" OFS "N"
                next
            }
            {
                key = \$association_predictor SUBSEP \$association_a1 SUBSEP \$association_a2
                if (key in seen_association) {
                    print "Duplicate Predictor+A1+A2 key in LDAK-KVIK association table: " \$association_predictor, \$association_a1, \$association_a2 > "/dev/stderr"
                    exit 1
                }
                seen_association[key] = 1
                if (!(key in sample_count)) {
                    print "No LDAK-KVIK summary row for Predictor+A1+A2 key: " \$association_predictor, \$association_a1, \$association_a2 > "/dev/stderr"
                    exit 1
                }
                first = 1
                for (i = 1; i <= NF; i++) {
                    if (header[i] != "Wald_Stat" && header[i] != "MAF") {
                        printf "%s%s", (first ? "" : OFS), \$i
                        first = 0
                    }
                }
                print OFS effect_allele_frequency[key] OFS sample_count[key]
            }
            END {
                if (!association_predictor) {
                    exit 1
                }
                for (key in sample_count) {
                    if (!(key in seen_association)) {
                        print "LDAK-KVIK summary key has no association row" > "/dev/stderr"
                        exit 1
                    }
                }
            }
        ' "\${summary}" "\${assoc}" > "\${harmonisation_input}"
    done
    """

    stub:
    prefix = task.ext.prefix ?: meta2.id
    """
    touch "${prefix}.step2.assoc"
    touch "${prefix}.step2.harmonisation.tsv"
    touch "${prefix}.step2.summaries"
    touch "${prefix}.step2.pvalues"
    touch "${prefix}.step2.log"
    """
}
