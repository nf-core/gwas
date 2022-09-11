missingness = [0.01,0.03,0.05]  // this is used by one of the templates

process ANALYZE_X {
     input:
       path(xchr)

     output:
        path(out), emit: x_analy_res_ch // batchReport

     script:
        def x = xchr[0].baseName
        def out = "x.pkl"
        //def f_hi_female = "" FIX
        //def f_lo_male "" FIX
    """
    xCheck.py \\
    --x ${x} \\
    --f_hi_female ${f_hi_female} \\
    --f_lo_male ${f_lo_male} \\
    --out ${out} 
    """
   }

