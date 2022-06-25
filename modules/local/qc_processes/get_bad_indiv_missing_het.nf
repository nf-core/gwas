// Find those who have bad heterozygosity
process get_bad_indivs_missing_het {
  memory other_mem_req
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.txt"

  input:
    tuple path(het), path(imiss)

  output:
    path(outfname), emit: failed_miss_het

  script:
    def base = het.baseName
    def outfname = "${base}-fail_het".replace(".","_")+".txt"

    template "select_miss_het_qcplink.py"
}
