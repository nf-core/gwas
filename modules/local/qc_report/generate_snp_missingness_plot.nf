process GENERATE_SNP_MISSINGNESS_PLOT {
  debug true

  input:
      path(lmissf)

  output:
      path("${output}"), emit: report_snpmiss_ch


  script:
    def input  = lmissf
    def base   = lmissf.baseName
    def label  = "SNPs"
    output = "${base}-snpmiss_plot".replace(".","_")+".pdf"
    """
    missPlot.py --input $lmissf --label "SNPs" --output $output
    """
}
