process ExtractRsIDChro{
    memory params.mem_req 
    input :
     set val(chro), file(gwas), file(rsinfo) from gwas_format_chro_rs
    output :
      set val(chro), file(gwas),file(outrs) into rsinfo_chro
    script :
      outrs="info_rs_"+chro+".rs"
    """
    zcat $rsinfo | extractrsid_bypos.py --file_chrbp $gwas --out_file $outrs --ref_file stdin --chr $chro --chro_ps ${params.poshead_chro_inforef} --bp_ps ${params.poshead_bp_inforef} --rs_ps ${params.poshead_rs_inforef} --a1_ps ${params.poshead_a1_inforef}  --a2_ps ${params.poshead_a2_inforef} 
    """ 
}
