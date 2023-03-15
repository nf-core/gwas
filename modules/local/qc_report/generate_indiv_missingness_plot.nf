process GENERATE_INDIV_MISSINGNESS_PLOT {
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
        path(imissf)

  output:
        path("${output}"), emit: report_indmisspdf_ch

  script:
        def input  = imissf
        def base   = imissf.baseName
        output = "${base}-indmiss_plot".replace(".","_")+".pdf"
         """
         missPlot.py --input $imissf --label samples --output $output
         """
}
