// Find those SNPs that have diff missingness in cases & controls
process findSnpExtremeDifferentialMissingness {
  memory other_mem_req

    input:
    path(clean_missing)
  echo true

    output:
    tuple val(base), path(failed), emit: bad_snps_ch
    path(failed), emit: report_diffmiss_ch
    path(failed), emit: skewsnps_ch
  script:
    cut_diff_miss=params.cut_diff_miss
    missing = clean_missing
    base     = missing.baseName.replace("-.*","").replace(".","_")
    probcol = 'EMP2'  // need to change if we don't use mperm
    failed   = "${base}-failed_diffmiss.snps"
    template "select_diffmiss_qcplink.py"
}
