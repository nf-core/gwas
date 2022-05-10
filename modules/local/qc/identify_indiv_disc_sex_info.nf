/* Process to identify individual discordant sex information.
 * results are put in the output directory
 * Also does HWE
 */
process identifyIndivDiscSexinfo {
  memory plink_mem_req
  input:
     path(plinks)

  publishDir params.output_dir, overwrite:true, mode:'copy'

  output:
    path(logfile),emit:  (report_failed_sex_ch, failed_sex_ch1)
    tuple path(imiss), path(lmiss),path(sexcheck_report),emit: batchrep_missing_ch
    path("${base}.hwe"), emit: hwe_stats_ch
  validExitStatus 0, 1
  script:
    base = plinks[0].baseName
    logfile= "${base}.badsex"
    sexcheck_report = "${base}.sexcheck"
    imiss  = "${base}.imiss"
    lmiss  = "${base}.lmiss"
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
