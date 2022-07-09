process REMOVE_QC_INDIVS {

  input:
        path(f_miss_het)
        path(rel_indivs)
        path (f_sex_check_f)
        path (poorgc)
        tuple path(bed), path(bim), path(fam)

  output:
        path("${out}.{bed,bim,fam}"), emit: qc3_ch

  script:
        def base = bed.baseName
        def out  = "${base}-c".replace(".","_")
            """
            cat $f_sex_check_f $rel_indivs $f_miss_het $poorgc | sort -k1 | uniq > failed_inds
            plink $K --bfile $base $sexinfo --remove failed_inds --make-bed --out $out
            mv failed_inds ${out}.irem
            """
}
