// Find differential missingness between cases and controls; also compute HWE scores
process CALCULATE_SNP_SKEW_STATUS {
  input:
        path(plinks)

  output:
        path("${base}.missing"), emit: clean_diff_miss_plot_ch
        path(mperm), emit: clean_diff_miss_ch
        path("${base}.hwe"), emit: hwe_scores_ch

  script:
        def base  = plinks[0].baseName
        def out   = base.replace(".","_")
        def mperm = "${base}.missing.mperm"
        def phe   = plinks[3]

        """
            cp $phe cc.phe
            plink --threads ${max_plink_cores} --autosome --bfile $base $sexinfo $diffpheno --test-missing mperm=10000 --hardy --out $out
            if ! [ -e $mperm ]; then
            echo "$mperm_header" > $mperm
            fi
        """
}
