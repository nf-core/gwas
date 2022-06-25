missingness = [0.01,0.03,0.05]  // this is used by one of the templates

if (extrasexinfo == "--must-have-sex") {

   /* Detailed analysis of X-chromosome */
   process get_x {

     memory other_mem_req

     input:
       path(plink)

      output:
        path("X*"),emit: X_chr_ch

      script:
        def base = plink[0].baseName

            """
            if [[ `grep  "^23" *bim`  ]];  then
            plink --bfile $base --chr 23 --geno 0.04 --make-bed --out X
            else
            echo ""
            echo "----------------------------------------"
            echo "There are no X-chromosome SNPs in this data "
            echo "it does not make sense to check for sex"
            echo "set sexinfo_available to false"
            echo "----------------------------------------"
            echo ""
            exit 23
            touch X.bed X.bim X.fam EMPTYX
            fi
            """
   }



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
} else {


  x_analy_res_ch = Channel.fromPath("0")


}
