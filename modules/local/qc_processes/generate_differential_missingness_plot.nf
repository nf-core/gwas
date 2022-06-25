process generate_differential_missingness_plot {
   memory other_mem_req

   input:
     file clean_missing

   output:
      file output

   script:
       def input = clean_missing
       def base  = clean_missing.baseName.replace(".","_").replace("-nd","")
       def output= "${base}-diff-snpmiss_plot.pdf"

       template "diffMiss.py"

 }
