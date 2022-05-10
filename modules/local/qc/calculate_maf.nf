process calculateMaf {
  memory plink_mem_req
  input:
    tuple  path(bed), path(bim), path(fam), path(log) from qc4C_ch

  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.frq"

  output:
    path "${base}.frq" into maf_plot_ch

  script:
    base = bed.baseName
    out  = base.replace(".","_")
    """
      plink --bfile $base $sexinfo  --freq --out $out
    """
}
