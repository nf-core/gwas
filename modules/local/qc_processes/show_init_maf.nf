process SHOW_INIT_MAF {
  input:
        path(freq)

  output:
        tuple path("${base}.pdf"), path("${base}.tex"), emit: report_initmaf_ch

  script:
        def base = freq.baseName+"-initmaf"
        def base = base.replace(".","_")

        template "showmaf.py"
}
