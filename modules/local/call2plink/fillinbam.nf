process fillInBim {  //  Deals with monomorphic or non-called SNPs
    input:
     file(inbim)
     file(strand)
     file(manifest)
    output:
     file("raw.bim"), emit: filled_bim_ch
    script:
       "fill_in_bim.py ${params.output_align} $strand $manifest $inbim  raw.bim"
  }
