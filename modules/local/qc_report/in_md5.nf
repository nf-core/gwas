process IN_MD5 {
  input:
     path(plink)
  output:
     path(out), emit: report_inpmd5_ch
  echo true
  script:
       pat = plink[0]
       bed = plink[1]
       bim = plink[2]
       fam = plink[3]
       out  = "${plink[0].baseName}.md5"
       template "md5.py"
}
