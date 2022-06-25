process calculate_sample_heterozygosity {
   memory plink_mem_req
   publishDir params.output_dir, overwrite:true, mode:'copy'

   input:
      path(nodups)


   output:
    tuple path("${hetf}.het"), path("${hetf}.imiss"), emit: hetero_check_ch
    path("${hetf}.imiss"), emit: missing_stats_ch

   script:
      def base = nodups[0].baseName
      def hetf = "${base}".replace(".","_")

    """
        plink --bfile $base  $sexinfo --het --missing  --out $hetf
    """
}
