
  process SAMPLE_SHEET {
    input:
     file(sheet)
    output:
     file("poorgc10.lst"), emit: poorgc10_ch
     file("plates"), emit: report_poorgc10_ch
    script:
     """
       mkdir -p plates
       sampleqc.py $sheet ${params.gc10} "${idpat}"  poorgc10.lst plates/crgc10.tex
      """
    }
