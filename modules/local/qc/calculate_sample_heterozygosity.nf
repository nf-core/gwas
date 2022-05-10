process calculateSampleHeterozygosity {
   memory plink_mem_req
   input:
      path(nodups) from qc2C_ch

   publishDir params.output_dir, overwrite:true, mode:'copy'
   output:
      tuple path("${hetf}.het"), path("${hetf}.imiss") into (hetero_check_ch, plot1_het_ch)
      path("${hetf}.imiss") into missing_stats_ch
   script:
      base = nodups[0].baseName
      hetf = "${base}".replace(".","_")
   """
     plink --bfile $base  $sexinfo --het --missing  --out $hetf
   """
}
