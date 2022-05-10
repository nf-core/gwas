
// Find those who have bad heterozygosity
process getBadIndivsMissingHet {
  memory other_mem_req
  input:
    tuple path(het), path(imiss)
  output:
    path(outfname), emit: failed_miss_het
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.txt"
  script:
    base = het.baseName
    outfname = "${base}-fail_het".replace(".","_")+".txt"
    template "select_miss_het_qcplink.py"
}
