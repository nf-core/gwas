process generateHwePlot {
  memory other_mem_req
  input:
    file unaff from unaff_hwe
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    file output into report_hwepdf_ch

  script:
    input  = unaff
    base   = unaff.baseName.replace(".","_")
    output = "${base}-hwe_plot.pdf"
    template "hweplot.py"
}
