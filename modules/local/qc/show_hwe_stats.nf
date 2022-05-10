process showHWEStats {
  memory other_mem_req
  input:
     path(hwe) from hwe_stats_ch
  output:
     tuple path("${base}.pdf"), path("${base}-qq.pdf"), path("${base}.tex") into report_inithwe_ch
  script:
    base = hwe.baseName+"-inithwe"
    base = base.replace(".","_")
    template "showhwe.py"
}
