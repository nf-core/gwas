process GENERATE_HWE_PLOT {
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
        path(unaff)

  output:
        path(output), emit: report_hwepdf_ch

  script:
        def input  = unaff
        def base   = unaff.baseName.replace(".","_")
        def output = "${base}-hwe_plot.pdf"

        template "hweplot.py"
}
