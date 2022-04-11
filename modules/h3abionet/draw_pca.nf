process drawPCA {
    memory other_mem_req
    input:
      tuple file(eigvals), file(eigvecs) from pcares
      file cc from cc2_ch
    output:
      tuple  file ("eigenvalue.pdf"), file(output) into report_pca_ch
    publishDir params.output_dir, overwrite:true, mode:'copy',pattern: "*.pdf"
    script:
      base=eigvals.baseName
      cc_fname = params.case_control
      // also relies on "col" defined above
      output="${base}-pca".replace(".","_")+".pdf"
      template "drawPCA.py"

}
