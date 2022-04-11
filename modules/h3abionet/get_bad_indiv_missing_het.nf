
// Find those who have bad heterozygosity
process getBadIndivsMissingHet {
  memory other_mem_req
  input:
    tuple file(het), file(imiss) from hetero_check_ch
  output:
    file(outfname) into (failed_miss_het, report_misshetremf_ch)
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.txt"
  script:
    base = het.baseName
    outfname = "${base}-fail_het".replace(".","_")+".txt"
    template "select_miss_het_qcplink.py"
}
