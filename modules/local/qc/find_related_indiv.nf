// run script to find a tuple of individuals we can remove to ensure no relatedness
//  Future - perhaps replaced with Primus
process findRelatedIndiv {
  memory other_mem_req
  input:
     path (missing) from ind_miss_ch2
     path (ibd_genome) from find_rel_ch
  output:
     path(outfname) into (related_indivs_ch1,related_indivs_ch2, report_related_ch)
  publishDir params.output_dir, overwrite:true, mode:'copy'
  script:
     base = missing.baseName
     outfname = "${base}-fail_IBD".replace(".","_")+".txt"
     template "removeRelInds.py"
}
