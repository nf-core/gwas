// Find those SNPs that have diff missingness in cases & controls
process FIND_SNP_EXTREME_DIFFERENTIAL_MISSINGNESS {
  debug true

  input:
        path(clean_missing)

  output:
        tuple val(base), path(failed), emit: bad_snps_ch
        path(failed), emit: skewsnps_ch

  script:
        def cut_diff_miss=params.cut_diff_miss
        def missing = clean_missing
        def base     = missing.baseName.replace("-.*","").replace(".","_")
        def probcol = 'EMP2'  // need to change if we don't use mperm
        def failed   = "${base}-failed_diffmiss.snps"
        """
        select_diffmiss_qcplink.py \\
        --missing ${missing}\\
        --probcol ${probcol}\\
        --cut_diff_miss ${cut_diff_miss}\\
        --failed ${failed}
        """
}
