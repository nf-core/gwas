process showInitMAF {
  memory other_mem_req
  input:
     file(freq) from init_freq_ch
  output:
     tuple file("${base}.pdf"), file("${base}.tex") into report_initmaf_ch
  script:
    base = freq.baseName+"-initmaf"
    base = base.replace(".","_")
    template "showmaf.py"
}
