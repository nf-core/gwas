repnames = ["dups","cleaned","misshet","mafpdf","snpmiss","indmisspdf","failedsex","misshetremf","diffmissP","diffmiss","pca","hwepdf","related","inpmd5","outmd5","batch"]



process produceReports {
  memory other_mem_req
  label 'latex'
  input:
    tuple file(orig), file (dupf) from report_dups_ch
    tuple file(cbed), file(cbim), file(cfam), file(ilog) from report_cleaned_ch
    file(missingvhetpdf) from report_misshet_ch
    file(mafpdf)         from report_mafpdf_ch
    file(snpmisspdf)     from report_snpmiss_ch
    file(indmisspdf)     from report_indmisspdf_ch
    file(fsex)           from report_failed_sex_ch
    file(misshetremf)    from report_misshetremf_ch
    file(diffmisspdf)    from report_diffmissP_ch
    file(diffmiss)       from report_diffmiss_ch
    tuple file(eigenvalpdf),file(pcapdf)         from report_pca_ch
    file(hwepdf)         from report_hwepdf_ch
    file(rel_indivs)     from report_related_ch
    file(inpmd5)         from report_inpmd5_ch
    file(outmd5)         from report_outmd5_ch
    tuple file(initmafpdf), file(initmaftex) from report_initmaf_ch
    tuple file(inithwepdf), file(inithweqqpdf), file(inithwetex) from report_inithwe_ch
    tuple file(qc1), file(irem)  from report_qc1_ch
    file(batch_tex)  from report_batch_report_ch
    file(poorgc)     from report_poorgc10_ch
    tuple file(bpdfs), file(bcsvs) from report_batch_aux_ch
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    file("${base}.pdf") into final_ch
   script:
     base = params.output
     config_text = getConfig()
     template "qcreport.py"
}
