process convertInVcf {
   memory plink_mem_req
   cpus max_plink_cores
   input :
     tuple file(bed), file(bim), file(fam), file (log) from qc4A_ch
   publishDir params.output_dir, overwrite:true, mode:'copy'
   output :
    file("${base}.vcf")
   script:
     base= bed.baseName
     """
     plink --bfile ${base} --threads ${max_plink_cores} --recode vcf --out $base
     """
}
