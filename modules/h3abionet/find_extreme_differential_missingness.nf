// Find those SNPs that have diff missingness in cases & controls
process findSnpExtremeDifferentialMissingness {
  memory other_mem_req
  input:
    file clean_missing from clean_diff_miss_ch2
  echo true
  output:
     tuple val(base), file(failed) into bad_snps_ch
     file(failed) into report_diffmiss_ch
     file(failed) into skewsnps_ch
  script:
    cut_diff_miss=params.cut_diff_miss
    missing = clean_missing
    base     = missing.baseName.replace("-.*","").replace(".","_")
    probcol = 'EMP2'  // need to change if we don't use mperm
    failed   = "${base}-failed_diffmiss.snps"
    template "select_diffmiss_qcplink.py"
}
