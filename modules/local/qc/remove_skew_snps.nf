process removeSkewSnps {
  memory plink_mem_req
  input:
    path (plinks)
    path(failed)
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    tuple path("${output}.bed"), path("${output}.bim"), path("${output}.fam"), path("${output}.log"), emit: qc4A_ch
  script:
  base = plinks[0].baseName
  output = params.output.replace(".","_")
  """
  plink $K --bfile $base $sexinfo --exclude $failed --make-bed --out $output
  """
}
