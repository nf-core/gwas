process MergeRsGwasChro{
    memory params.mem_req 
    input :
      set val(chro), file(gwas),file(chrors), file(bed), file(bim), file(fam) from rsinfo_chroplk
    output :
      file(outmerge) into gwas_rsmerge
    script :
     outmerge="merge_"+chro+".gwas"
     bfileopt= (params.input_pat!="" || params.input_dir!="") ?  " --bfile "+bed.baseName : ""
     Nheadopt=(params.head_N!="") ? " --N_head ${params.head_N} " : ""
     Freqheadopt=(params.head_freq!="") ? " --freq_head ${params.head_freq} " : ""

     NheadNewopt=(params.headnew_N!="") ? " --Nnew_head ${params.headnew_N} " : ""
     FreqNewheadopt=(params.headnew_freq!="") ? " --freqnew_head ${params.headnew_freq} " : ""
     """
     mergeforrs.py --input_gwas $gwas --input_rs $chrors  --out_file $outmerge --chro_head  ${params.headnew_chr} --bp_head  ${params.headnew_bp} --rs_head ${params.headnew_rs} --chro $chro $bfileopt  $Nheadopt $Freqheadopt $NheadNewopt $FreqNewheadopt  --a1_head ${params.headnew_A1} --a2_head  ${params.headnew_A2}
     """

}


