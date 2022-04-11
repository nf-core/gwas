process removeQCIndivs {
  memory plink_mem_req
  input:
    file(f_miss_het)     from failed_miss_het
    file(rel_indivs)     from related_indivs_ch1
    file (f_sex_check_f) from failed_sex_ch1
    file (poorgc)        from poorgc10_ch
    tuple file(bed), file(bim), file(fam) from qc2D_ch
  output:
     file("${out}.{bed,bim,fam}") into\
        (qc3A_ch, qc3B_ch)
  script:
   base = bed.baseName
   out  = "${base}-c".replace(".","_")
    """
     cat $f_sex_check_f $rel_indivs $f_miss_het $poorgc | sort -k1 | uniq > failed_inds
     plink $K --bfile $base $sexinfo --remove failed_inds --make-bed --out $out
     mv failed_inds ${out}.irem
  """
}
