// Combine  all plinks and possibly remove errprs
 process combineBed {
   input:
     file(bed) 
     file(bim) 
     file(fam)
   output:
     tuple file("raw.bed"), file("raw.fam"), file("raw.log"), emit: plink_src
     file("rawraw.bim"), emit: fill_in_bim_ch
   script:
    """
    ls *.bed | sort > beds
    ls *.bim | sort >  bims
    ls *.fam | sort >  fams
    paste beds bims fams >  mergelist
    plink --merge-list mergelist  --make-bed --out raw
    mv raw.bim rawraw.bim
    """
 } 
