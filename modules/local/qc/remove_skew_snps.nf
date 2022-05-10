process removeSkewSnps {
  memory plink_mem_req
  input:
    path (plinks) from qc3B_ch
    path(failed) from skewsnps_ch
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    tuple path("${output}.bed"), path("${output}.bim"), path("${output}.fam"), path("${output}.log") \
      into (qc4A_ch, qc4B_ch, qc4C_ch,  report_cleaned_ch)
  script:
  base = plinks[0].baseName
  output = params.output.replace(".","_")
  """
  plink $K --bfile $base $sexinfo --exclude $failed --make-bed --out $output
  """
}
