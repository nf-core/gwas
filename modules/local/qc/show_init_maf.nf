process showInitMAF {
  memory other_mem_req
  input:
     path(freq) from init_freq_ch
  output:
     tuple path("${base}.pdf"), path("${base}.tex") into report_initmaf_ch
  script:
    base = freq.baseName+"-initmaf"
    base = base.replace(".","_")
    template "showmaf.py"
}
