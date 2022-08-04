process GET_X {

  input:
        path(plink)

  output:
        path("X*"), emit: x_chr_ch

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
