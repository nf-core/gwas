process showInitMAF {
  memory other_mem_req
  input:
     path(freq)
  output:
    tuple path("${base}.pdf"), path("${base}.tex"), emit: report_initmaf_ch
  script:
    base = freq.baseName+"-initmaf"
    base = base.replace(".","_")
    template "showmaf.py"
}
