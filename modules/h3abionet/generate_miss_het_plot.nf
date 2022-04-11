process generateMissHetPlot {
  memory other_mem_req
  input:
    tuple file(het), file(imiss) from plot1_het_ch
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    file(output) into report_misshet_ch
  script:
    base = imiss.baseName
    output  = "${base}-imiss-vs-het".replace(".","_")+".pdf"
    template "missHetPlot.py"
}
