process GENERATE_MAF_PLOT {

  input:
        path(input)

  output:
        path(output)

  script:
        def base    = input.baseName
        def output  = "${base}-maf_plot.pdf"

        template "mafplot.py"
}
