// Generate MD5 sums of output files
process outMD5 {
  memory other_mem_req
  input:
     tuple file(bed), file(bim), file(fam), file(log) from qc4B_ch
  output:
     file(out) into report_outmd5_ch
  echo true
  script:
       out  = "${bed.baseName}.md5"
       template "md5.py"
}
