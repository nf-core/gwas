process SHOW_HWE_STATS {

  input:
        path(hwe)

  output:
        tuple path("${base}.pdf"), path("${base}-qq.pdf"), path("${base}.tex"), emit: report_inithwe_ch

  script:
        def base = hwe.baseName+"-inithwe"
        def base = base.replace(".","_")

        template "showhwe.py"
}
