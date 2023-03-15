process DRAW_PCA {
    publishDir params.output_dir, overwrite:true, mode:'copy',pattern: "*.pdf"

    input:
      tuple file(eigvals), file(eigvecs)
      file(cc)

    output:
      tuple  file ("eigenvalue.pdf"), file(output), emit: report_pca_ch

    script:
      def base=eigvals.baseName
      def cc_fname = params.case_control
      def col= params.case_control_col
      def output="${base}-pca".replace(".","_")+".pdf"
      """
      drawPCA.py --base $base \\
      --cc $cc \\
      --cc_fname $cc_fname \\
      --col $col \\
      --eigvals $eigvals \\
      --eigvecs $eigvecs \\
      --output $output

      """

}
