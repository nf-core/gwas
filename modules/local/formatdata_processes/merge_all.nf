process MergeAll{
   memory params.mem_req 
   input :
      file(allfile) from gwas_rsmerge_all
   publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
   output :
      file(fileout) 
   script :
     file1=allfile[0]
     listefiles=allfile.join(" ")
     fileout=params.output
     """
     head -1 $file1 > $fileout
     ls $listefiles  | xargs -n 1 tail -n +2 >> $fileout
     """

}



