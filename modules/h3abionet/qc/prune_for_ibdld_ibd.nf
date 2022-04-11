if (params.high_ld_regions_fname != "") {

  ldreg_ch=Channel.fromPath(params.high_ld_regions_fname)

  process pruneForIBDLD {
    cpus max_plink_cores
    memory plink_mem_req
    input:
      file plinks from qc2B_ch
      file ldreg  from ldreg_ch
      publishDir params.output_dir, overwrite:true, mode:'copy'
    output:
      file "${outf}.genome" into (find_rel_ch,batch_rel_ch)
    script:
      base   = plinks[0].baseName
      outf   =  base.replace(".","_")
      range = " --exclude range $ldreg"
      """
       plink --bfile $base --threads $max_plink_cores --autosome $sexinfo $range --indep-pairwise 60 5 0.2 --out ibd
       plink --bfile $base --threads $max_plink_cores --autosome $sexinfo --extract ibd.prune.in --genome --out ibd_prune
       plink --bfile $base --threads $max_plink_cores --autosome $sexinfo --extract ibd.prune.in --genome --min $pi_hat --out $outf
       echo LD
       """
  }
} else {

  process pruneForIBD {
    cpus max_plink_cores
    memory plink_mem_req
    input:
      file plinks from qc2B_ch
      publishDir params.output_dir, overwrite:true, mode:'copy'
    output:
      file "${outf}.genome" into (find_rel_ch,batch_rel_ch)
    script:
      base   = plinks[0].baseName
      outf   =  base.replace(".","_")
      """
       plink --bfile $base --threads $max_plink_cores --autosome $sexinfo  --indep-pairwise 60 5 0.2 --out ibd
       #plink --bfile $base --threads $max_plink_cores --autosome $sexinfo --extract ibd.prune.in --genome --out ibd_prune
       plink --bfile $base --threads $max_plink_cores --autosome $sexinfo --extract ibd.prune.in --genome --min $pi_hat --out $outf
       echo NO_LD
       """
  }
}