process generate_hwe_plot {
  memory other_mem_req
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
    file unaff from unaff_hwe

  output:
    file output into report_hwepdf_ch

  script:
    def input  = unaff
    def base   = unaff.baseName.replace(".","_")
    def output = "${base}-hwe_plot.pdf"

    template "hweplot.py"
}
