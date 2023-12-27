process ExtractChroGWAS{
    memory params.mem_req 
    input :
      file(gwas) from gwas_chrolist_ext
    each chro from  chrolist2
    output :
      set val(chro), file(gwas_out) into gwas_format_chro
    script :
      gwas_out=gwas.baseName+"_"+chro+".gwas"
      sep=(params.sep!="") ?  "" : ""
      infofile="Chro:${params.head_chr}:${params.headnew_chr},Pos:${params.head_bp}:${params.headnew_bp},A2:${params.head_A2}:${params.headnew_A2},A1:${params.head_A1}:${params.headnew_A1},af:${params.head_freq}:${params.headnew_freq},Beta:${params.head_beta}:${params.headnew_beta},Se:${params.head_se}:${params.headnew_se},Pval:${params.head_pval}:${params.headnew_pval},N:${params.head_N}:${params.headnew_N}"
      """
      extractandformat_gwas.py --input_file $gwas --out_file ${gwas_out} --chr $chro --info_file $infofile
      """
}


