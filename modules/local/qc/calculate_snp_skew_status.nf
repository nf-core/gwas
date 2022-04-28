mperm_header=" CHR                               SNP         EMP1         EMP2 "

// Find differential missingness between cases and controls; also compute HWE scores
process calculateSnpSkewStatus {
  memory plink_mem_req
  cpus max_plink_cores
  input:
    file(plinks) from qc3A_ch.combine(cc_ch)
  output:
    file "${base}.missing" into clean_diff_miss_plot_ch1
    file mperm into clean_diff_miss_ch2
    file "${base}.hwe" into hwe_scores_ch
  script:
   base  = plinks[0].baseName
   out   = base.replace(".","_")
   mperm = "${base}.missing.mperm"
   phe   = plinks[3]
   """
    cp $phe cc.phe
    plink --threads ${max_plink_cores} --autosome --bfile $base $sexinfo $diffpheno --test-missing mperm=10000 --hardy --out $out
    if ! [ -e $mperm ]; then
       echo "$mperm_header" > $mperm
    fi
   """
}
