// We do PCA on qc2 data because relatively few SNPs and individuals will be removed later and
// this is an expensive operation so we start early. Similarly for computing relatedness
process compPCA {
   cpus max_plink_cores
   memory plink_mem_req
   input:
      file plinks from qc2A_ch
   output:
      tuple file ("${prune}.eigenval"), file("${prune}.eigenvec") into (pcares, pcares1)
      tuple file ("${prune}.bed"), file("${prune}.bim"), file("${prune}.fam") into out_only_pcs_ch
   publishDir "${params.output_dir}/pca", overwrite:true, mode:'copy',pattern: "${prune}*"
   script:
      base = plinks[0].baseName
      prune= "${base}-prune".replace(".","_")
     """
     plink --bfile ${base} --indep-pairwise 100 20 0.2 --out check
     plink --bfile ${base} --extract check.prune.in --make-bed --out $prune
     /bin/rm check*
     plink --bfile ${prune} --pca --out $prune
     """
}
