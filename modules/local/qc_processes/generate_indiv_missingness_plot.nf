process GENERATE_INDIV_MISSINGNESS_PLOT {
  memory other_mem_req
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
      path(imissf)

  output:
      file(output), emit: report_indmisspdf_ch

  script:
    def input  = imissf
    def base   = imissf.baseName
    def label  = "samples"
    def output = "${base}-indmiss_plot".replace(".","_")+".pdf"

    template "missPlot.py"
}
