process GENERATE_MISS_HET_PLOT {
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
        tuple file(het), file(imiss)

  output:
        path(output), emit: report_misshet_ch

  script:
        def base = imiss.baseName
        def output  = "${base}-imiss-vs-het".replace(".","_")+".pdf"

        template "missHetPlot.py"
}
