// We do PCA on qc2 data because relatively few SNPs and individuals will be removed later and
// this is an expensive operation so we start early. Similarly for computing relatedness
process COMP_PCA {
   publishDir "${params.output_dir}/pca", overwrite:true, mode:'copy',pattern: "${prune}*"

   input:
        path(plinks)

   output:
        tuple path("${prune}.eigenval"), path("${prune}.eigenvec"), emit: pcares
        tuple path("${prune}.bed"), path("${prune}.bim"), path("${prune}.fam"), emit: out_only_pcs_ch

   script:
        def base = plinks[0].baseName
        def prune = ("${base}-prune").replace(".","_")

        """
        plink --bfile ${base} --indep-pairwise 100 20 0.2 --out check
        plink --bfile ${base} --extract check.prune.in --make-bed --out $prune
        /bin/rm check*
        plink --bfile ${prune} --pca --out $prune
        """
}
