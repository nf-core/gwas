process generateIndivMissingnessPlot {
  memory other_mem_req
  input:
      file(imissf) from ind_miss_ch1
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    file(output) into report_indmisspdf_ch
  script:
    input  = imissf
    base   = imissf.baseName
    label  = "samples"
    output = "${base}-indmiss_plot".replace(".","_")+".pdf"
    template "missPlot.py"
}
