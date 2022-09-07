process GET_INIT_MAF {

  input:
        path(plink)

  output:
        path("${newbase}.frq"), emit: init_freq_ch

  script:
        def base = plink[0].baseName
        def newbase = base.replace(".","_")

        """
        plink --bfile $base --freq --out $newbase
        """
}