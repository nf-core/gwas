process removeQCIndivs {
  memory plink_mem_req
  input:
    path(f_miss_het)     from failed_miss_het
    path(rel_indivs)     from related_indivs_ch1
    path (f_sex_check_f) from failed_sex_ch1
    path (poorgc)        from poorgc10_ch
    tuple path(bed), path(bim), path(fam) from qc2D_ch
  output:
     path("${out}.{bed,bim,fam}") into\
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
