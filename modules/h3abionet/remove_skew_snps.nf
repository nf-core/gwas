process removeSkewSnps {
  memory plink_mem_req
  input:
    file (plinks) from qc3B_ch
    file(failed) from skewsnps_ch
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    tuple file("${output}.bed"), file("${output}.bim"), file("${output}.fam"), file("${output}.log") \
      into (qc4A_ch, qc4B_ch, qc4C_ch,  report_cleaned_ch)
  script:
  base = plinks[0].baseName
  output = params.output.replace(".","_")
  """
  plink $K --bfile $base $sexinfo --exclude $failed --make-bed --out $output
  """
}
