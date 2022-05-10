// Generate MD5 sums of output paths
process outMD5 {
  memory other_mem_req
  input:
     tuple path(bed), path(bim), path(fam), path(log)

  output:
    path(out), emit: report_outmd5_ch
  echo true
  script:
       out  = "${bed.baseName}.md5"
       template "md5.py"
}
