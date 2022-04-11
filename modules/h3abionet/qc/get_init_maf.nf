process getInitMAF {
  memory plink_mem_req
  input:
     file(plink) from qc1C_ch
  output:
     file("${newbase}.frq") into init_freq_ch
  script:
    base = plink[0].baseName
    newbase = base.replace(".","_")
    """
    plink --bfile $base --freq --out $newbase
    """
}
