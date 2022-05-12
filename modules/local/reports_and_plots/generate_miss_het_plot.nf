process generateMissHetPlot {
  memory other_mem_req
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
    tuple file(het), file(imiss)

  output:
    file(output)

  script:
    def base = imiss.baseName
    def output  = "${base}-imiss-vs-het".replace(".","_")+".pdf"

    template "missHetPlot.py"
}
