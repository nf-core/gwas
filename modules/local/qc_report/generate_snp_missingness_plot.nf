process GENERATE_SNP_MISSINGNESS_PLOT {
  debug true

  input:
      path(lmissf)

  output:
      path(output), emit: report_snpmiss_ch


  script:
    def input  = lmissf
    def base   = lmissf.baseName
    def label  = "SNPs"
    def output = "${base}-snpmiss_plot".replace(".","_")+".pdf"

    template "missPlot.py"
}
