process generateSnpMissingnessPlot {
  memory other_mem_req
  input:
      file(lmissf) from snp_miss_ch
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
     file(output) into report_snpmiss_ch

  echo true
  script:
    input  = lmissf
    base   = lmissf.baseName
    label  = "SNPs"
    output = "${base}-snpmiss_plot".replace(".","_")+".pdf"
    template "missPlot.py"
}
