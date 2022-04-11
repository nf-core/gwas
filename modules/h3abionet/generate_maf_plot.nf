process generateMafPlot {
  memory other_mem_req
  input:
    file input from maf_plot_ch
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    file(output) into report_mafpdf_ch

  script:
    base    = input.baseName
    output  = "${base}-maf_plot.pdf"
    template "mafplot.py"
}
