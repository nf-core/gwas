process REMOVE_SKEW_SNPS {

  input:
        path (plinks)
        path(failed)

  output:
        tuple path("${output}.bed"), path("${output}.bim"), path("${output}.fam"), path("${output}.log"), emit: qc4_ch

  script:
        def base = plinks[0].baseName
        def output = params.output.replace(".","_")
        """
        plink $K --bfile $base $sexinfo --exclude $failed --make-bed --out $output
        """
}
