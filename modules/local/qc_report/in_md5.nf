process IN_MD5 {
  input:
     tuple val(pattern), path(bed), path(bim), path(fam)
  output:
     path(out), emit: report_inpmd5_ch
  echo true
  script:
       out  = "${pattern}.md5"
       """
       md5.py --bed ${bed} --bim ${bim} --fam ${fam} --out ${out}
       """
}
