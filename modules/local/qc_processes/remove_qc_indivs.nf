process REMOVE_QC_INDIVS {

  input:
        path(failed_miss_het)
        path(related_indivs)
        path(failed_sex)
        path(poor_gc10)
        tuple path(bed), path(bim), path(fam)

  output:
        path("${out}.{bed,bim,fam}"), emit: qc3_ch

  script:
        def base = bed.baseName
        def out  = "${base}-c".replace(".","_")
            """
            cat $failed_sex $related_indivs $failed_miss_het $poor_gc10 | sort -k1 | uniq > failed_inds
            plink $K --bfile $base $sexinfo --remove failed_inds --make-bed --out $out
            mv failed_inds ${out}.irem
            """
}
