process generateMafPlot {
  memory other_mem_req

  input:
    file input

  output:
    file(output)

  script:
    def base    = input.baseName
    def output  = "${base}-maf_plot.pdf"

    template "mafplot.py"
}
