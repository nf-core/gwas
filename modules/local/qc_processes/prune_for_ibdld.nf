process PRUNE_FOR_IBDLD {

  input:
        path plinks
        path ldreg

  output:
        path "${outf}.genome", emit: find_rel_ch

  script:
        def base   = plinks[0].baseName
        def outf   =  base.replace(".","_")
        def range = " --exclude range $ldreg"

        """
        plink --bfile $base --threads $max_plink_cores --autosome $sexinfo $range --indep-pairwise 60 5 0.2 --out ibd
        plink --bfile $base --threads $max_plink_cores --autosome $sexinfo --extract ibd.prune.in --genome --out ibd_prune
        plink --bfile $base --threads $max_plink_cores --autosome $sexinfo --extract ibd.prune.in --genome --min $pi_hat --out $outf
        echo LD
        """
}
