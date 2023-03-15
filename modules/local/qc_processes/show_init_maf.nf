process SHOW_INIT_MAF {
  input:
        path(freq)

  output:
        tuple path("${base}.pdf"), path("${base}.tex"), emit: report_initmaf_ch

  script:
        base = (freq.baseName+"-initmaf").replace(".","_")

        """
        showmaf.py \\
        --freq ${freq} \\
        --base ${base}
        """
}
