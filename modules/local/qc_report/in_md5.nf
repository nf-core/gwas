process IN_MD5 {
  input:
     path(plink)
  output:
     path(out), emit: report_inpmd5_ch
  echo true
  script:
       bed = plink[0]
       bim = plink[1]
       fam = plink[2]
       out  = "${plink[0].baseName}.md5"
       template "md5.py"
}
