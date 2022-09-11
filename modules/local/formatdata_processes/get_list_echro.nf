process getListeChro{
        input :
          file(gwas_res) from gwas_chrolist
        output :
          file("filechro") into chrolist
        script:
         sep=(params.sep!="") ?  " --sep ${params.sep}" : ""
         """
         extractlistchro.py --input_file $gwas_res --chro_header ${params.head_chr} $sep > filechro
        """
}


