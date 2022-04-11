process showHWEStats {
  memory other_mem_req
  input:
     file(hwe) from hwe_stats_ch
  output:
     tuple file("${base}.pdf"), file("${base}-qq.pdf"), file("${base}.tex") into report_inithwe_ch
  script:
    base = hwe.baseName+"-inithwe"
    base = base.replace(".","_")
    template "showhwe.py"
}
