missingness = [0.01,0.03,0.05]  // this is used by one of the templates

process ANALYSE_X {
     input:
       path(xchr)

     output:
        path(out), emit: x_analy_res_ch // batchReport

     script:
        def x = xchr[0].baseName
        def out = "x.pkl"
        template "xCheck.py"
   }

