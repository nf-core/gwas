process GENERATE_DIFFERENTIAL_MISSINGNESS_PLOT {

   input:
        path(clean_missing)

   output:
        path(output) //report_diffmissP_ch

   script:
       def input = clean_missing
       def base  = clean_missing.baseName.replace(".","_").replace("-nd","")
       def output= "${base}-diff-snpmiss_plot.pdf"

       template "diffMiss.py"

 }
