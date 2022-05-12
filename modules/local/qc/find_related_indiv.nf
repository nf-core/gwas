// run script to find a tuple of individuals we can remove to ensure no relatedness
//  Future - perhaps replaced with Primus
process findRelatedIndiv {
  memory other_mem_req
  publishDir params.output_dir, overwrite:true, mode:'copy'

  input:
     path (missing)
     path (ibd_genome)

  output:
     path(outfname), emit: related_indivs_ch1

  script:
     def base = missing.baseName
     def outfname = "${base}-fail_IBD".replace(".","_")+".txt"
     template "removeRelInds.py"
}
