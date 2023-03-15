process SHOW_HWE_STATS {

  input:
        path(hwe)

  output:
        tuple path("${base}.pdf"), path("${base}-qq.pdf"), path("${base}.tex"), emit: report_inithwe_ch

  script:
        base = (hwe.baseName+"-inithwe").replace(".","_")
        
        """
        showhwe.py --hwe $hwe --base $base
        """
}
