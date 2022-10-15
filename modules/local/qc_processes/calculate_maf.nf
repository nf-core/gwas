process CALCULATE_MAF {
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.frq"

  input:
    tuple  path(bed), path(bim), path(fam), path(log)


  output:
    path "${base}.frq", emit: maf_plot_ch

  script:
    def base = bed.baseName
    def out  = base.replace(".","_")
    """
      plink --bfile $base $sexinfo  --freq --out $out
    """
}
