process SAMPLESHEET {
  input:
        path(sheet)

  output:
        path("poorgc10.lst"), emit into poorgc10_ch
        path("plates"), emit: report_poorgc10_ch

  script:
        """
        mkdir -p plates
        sampleqc.py $sheet ${params.gc10} "${idpat}"  poorgc10.lst plates/crgc10.tex
        """
}
