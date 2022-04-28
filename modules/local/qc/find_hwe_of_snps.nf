
// Find HWE scores of each SNP
process findHWEofSNPs {
  memory other_mem_req
  input:
     file hwe from hwe_scores_ch
  output:
     file output  into unaff_hwe

  script:
    base   = hwe.baseName.replace(".","_")
    output = "${base}-unaff.hwe"
    """
      head -1 $hwe > $output
      grep 'UNAFF' $hwe >> $output
    """
}
