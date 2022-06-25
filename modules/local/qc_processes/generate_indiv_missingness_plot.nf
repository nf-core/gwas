process generate_indiv_missingness_plot {
  memory other_mem_req
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
      file(imissf) from ind_miss_ch1

  output:
    file(output) into report_indmisspdf_ch

  script:
    def input  = imissf
    def base   = imissf.baseName
    def label  = "samples"
    def output = "${base}-indmiss_plot".replace(".","_")+".pdf"

    template "missPlot.py"
}
