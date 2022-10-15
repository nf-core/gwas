process DRAW_PCA {
    publishDir params.output_dir, overwrite:true, mode:'copy',pattern: "*.pdf"

    input:
      tuple file(eigvals), file(eigvecs)
      file(cc)

    output:
      tuple  file ("eigenvalue.pdf"), file(output) into report_pca_ch

    script:
      def base=eigvals.baseName
      def cc_fname = params.case_control
      // also relies on "col" defined above
      def output="${base}-pca".replace(".","_")+".pdf"

      template "drawPCA.py"

}
