process generateSnpMissingnessPlot {
  memory other_mem_req
  echo true

  input:
      file(lmissf) from snp_miss_ch

  output:
     file(output) into report_snpmiss_ch


  script:
    def input  = lmissf
    def base   = lmissf.baseName
    def label  = "SNPs"
    def output = "${base}-snpmiss_plot".replace(".","_")+".pdf"

    template "missPlot.py"
}
