/* Process to identify individual discordant sex information.
 * results are put in the output directory
 * Also does HWE
 */
process identifyIndivDiscSexinfo {
  memory plink_mem_req
  validExitStatus 0, 1

  input:
     path(plinks)


  output:
    path(logfile),emit:  (report_failed_sex_ch, failed_sex_ch1)
    tuple path(imiss), path(lmiss),path(sexcheck_report),emit: batchrep_missing_ch
    path("${base}.hwe"), emit: hwe_stats_ch

  script:
    def base = plinks[0].baseName
    def logfile= "${base}.badsex"
    def sexcheck_report = "${base}.sexcheck"
    def imiss  = "${base}.imiss"
    def lmiss  = "${base}.lmiss"

    if (params.sexinfo_available == true)
    """
       plink $K --bfile $base --hardy --check-sex $f_hi_female $f_lo_male --missing  --out $base
       head -n 1 ${base}.sexcheck > $logfile
       grep  'PROBLEM' ${base}.sexcheck >> $logfile
    """
    else
     """
     plink --bfile $base  --hardy --missing  --out $base
     echo 'FID IID STATUS' > $sexcheck_report
     echo 'No sex information'  > $logfile
     """
}
