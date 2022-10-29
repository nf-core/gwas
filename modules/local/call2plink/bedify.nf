process bedify {
    
   input:
     set file(ped), file(map)
   output:
     file("${base}.bed"), emit: bed_ch
     file ("${base}.bim"), emit: bim_ch
     file("${base}.fam"), emit: fam_ch
   script:
      def base = ped.baseName
      """
      echo $base

      plink --file $base --no-fid --make-bed --out $base
      """
  }
