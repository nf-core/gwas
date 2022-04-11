process generateDifferentialMissingnessPlot {
   memory other_mem_req
   input:
     file clean_missing from clean_diff_miss_plot_ch1
   publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"
   output:
      file output into report_diffmissP_ch
   script:
       input = clean_missing
       base  = clean_missing.baseName.replace(".","_").replace("-nd","")
       output= "${base}-diff-snpmiss_plot.pdf"
       template "diffMiss.py"

 }
