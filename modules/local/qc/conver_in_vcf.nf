process convert_in_vcf {
   memory plink_mem_req
   cpus max_plink_cores
   publishDir params.output_dir, overwrite:true, mode:'copy'

   input :
     tuple path(bed), path(bim), path(fam), path (log)

   output :
    path("${base}.vcf")

   script:
     def base= bed.baseName
     """
     plink --bfile ${base} --threads ${max_plink_cores} --recode vcf --out $base
     """
}
