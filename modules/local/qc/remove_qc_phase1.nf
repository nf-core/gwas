process removeQCPhase1 {
  memory plink_mem_req
  input:
    tuple path(bed), path(bim), path(fam)
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    path("${output}*.{bed,bim,fam}"), emit: qc2A_ch
    tuple path("qc1.out"), path("${output}.irem"),emit: report_qc1_ch
  script:
     base=bed.baseName
     output = "${base}-c".replace(".","_")
     """
     # remove really realy bad SNPs and really bad individuals
     plink $K --autosome --bfile $base $sexinfo --mind 0.1 --geno 0.1 --make-bed --out temp1
     plink $K --bfile temp1  $sexinfo --mind $params.cut_mind --make-bed --out temp2
     /bin/rm temp1*
     plink $K --bfile temp2  $sexinfo --geno $params.cut_geno --make-bed --out temp3
     /bin/rm temp2*
     plink $K --bfile temp3  $sexinfo --maf $params.cut_maf --make-bed --out temp4
     /bin/rm temp3*
     plink $K --bfile temp4  $sexinfo --hwe $params.cut_hwe --make-bed  --out $output
     /bin/rm temp4*
     cat *log > logfile
     touch tmp.irem
     cat *.irem > ${output}.irem
     qc1logextract.py logfile ${output}.irem > qc1.out
  """
}
