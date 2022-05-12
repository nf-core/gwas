// Find HWE scores of each SNP
process findHWEofSNPs {
  memory other_mem_req

  input:
     path hwe
  output:
    path output, emit: unaff_hwe

  script:
    def base   = hwe.baseName.replace(".","_")
    def output = "${base}-unaff.hwe"

    """
      head -1 $hwe > $output
      grep 'UNAFF' $hwe >> $output
    """
}
