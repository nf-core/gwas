process getInitMAF {
  memory plink_mem_req
  input:
     path(plink)
  output:
    path("${newbase}.frq"), emit: init_freq_ch
  script:
    base = plink[0].baseName
    newbase = base.replace(".","_")
    """
    plink --bfile $base --freq --out $newbase
    """
}
