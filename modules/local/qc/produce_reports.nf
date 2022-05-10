repnames = ["dups","cleaned","misshet","mafpdf","snpmiss","indmisspdf","failedsex","misshetremf","diffmissP","diffmiss","pca","hwepdf","related","inpmd5","outmd5","batch"]



process produceReports {
  memory other_mem_req
  label 'latex'
  input:
    tuple path(orig), path (dupf)
    tuple path(cbed), path(cbim), path(cfam), path(ilog)
    path(missingvhetpdf)
    path(mafpdf)
    path(snpmisspdf)
    path(indmisspdf)
    path(fsex)
    path(misshetremf)
    path(diffmisspdf)
    path(diffmiss)
    tuple path(eigenvalpdf),path(pcapdf)
    path(hwepdf)
    path(rel_indivs)
    path(inpmd5)
    path(outmd5)
    tuple path(initmafpdf), path(initmaftex)
    tuple path(inithwepdf), path(inithweqqpdf), path(inithwetex)
    tuple path(qc1), path(irem)
    path(batch_tex)
    path(poorgc)
    tuple path(bpdfs), path(bcsvs)
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    path("${base}.pdf"), emit: final_ch
   script:
     base = params.output
     config_text = getConfig()
     template "qcreport.py"
}
