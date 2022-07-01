process noSampleSheet {
    output:
     file("poorgc10.lst"), emit: poorgc10_ch
     file("plates"), emit: report_poorgc10_ch
    script:
     """
       mkdir -p plates
       sampleqc.py 0 0 0 poorgc10.lst plates/crgc10.tex
      """
    }
