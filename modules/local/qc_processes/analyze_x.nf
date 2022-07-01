missingness = [0.01,0.03,0.05]  // this is used by one of the templates

   process analyse_x {
     memory other_mem_req
     input:
       path(xchr)
     output:
        path(out), emit: x_analy_res_ch // batchReport
     script:
	x = xchr[0].baseName
	def out = "x.pkl"
	template "xCheck.py"
   }

