repnames = ["dups","cleaned","misshet","mafpdf","snpmiss","indmisspdf","failedsex","misshetremf","diffmissP","diffmiss","pca","hwepdf","related","inpmd5","outmd5","batch"]



process produceReports {
  memory other_mem_req
  label 'latex'
  input:
    tuple path(orig), path (dupf) from report_dups_ch
    tuple path(cbed), path(cbim), path(cfam), path(ilog) from report_cleaned_ch
    path(missingvhetpdf) from report_misshet_ch
    path(mafpdf)         from report_mafpdf_ch
    path(snpmisspdf)     from report_snpmiss_ch
    path(indmisspdf)     from report_indmisspdf_ch
    path(fsex)           from report_failed_sex_ch
    path(misshetremf)    from report_misshetremf_ch
    path(diffmisspdf)    from report_diffmissP_ch
    path(diffmiss)       from report_diffmiss_ch
    tuple path(eigenvalpdf),path(pcapdf)         from report_pca_ch
    path(hwepdf)         from report_hwepdf_ch
    path(rel_indivs)     from report_related_ch
    path(inpmd5)         from report_inpmd5_ch
    path(outmd5)         from report_outmd5_ch
    tuple path(initmafpdf), path(initmaftex) from report_initmaf_ch
    tuple path(inithwepdf), path(inithweqqpdf), path(inithwetex) from report_inithwe_ch
    tuple path(qc1), path(irem)  from report_qc1_ch
    path(batch_tex)  from report_batch_report_ch
    path(poorgc)     from report_poorgc10_ch
    tuple path(bpdfs), path(bcsvs) from report_batch_aux_ch
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    path("${base}.pdf") into final_ch
   script:
     base = params.output
     config_text = getConfig()
     template "qcreport.py"
}
