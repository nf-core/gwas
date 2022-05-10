process convertInVcf {
   memory plink_mem_req
   cpus max_plink_cores
   input :
     tuple path(bed), path(bim), path(fam), path (log) from qc4A_ch
   publishDir params.output_dir, overwrite:true, mode:'copy'
   output :
    path("${base}.vcf")
   script:
     base= bed.baseName
     """
     plink --bfile ${base} --threads ${max_plink_cores} --recode vcf --out $base
     """
}
