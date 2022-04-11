process batchProc {
  memory plink_mem_req
  input:
    tuple file(eigenval), file(eigenvec) from pcares1
    tuple file(imiss), file(lmiss), file(sexcheck_report) from batchrep_missing_ch
    file "pheno.phe" from phenotype_ch    // staged input file
    file "batch.phe" from batch_ch        // staged input file
    file genome    from batch_rel_ch    // pruneForIBD
    file pkl       from x_analy_res_ch  // analyseX
    file rem_indivs from related_indivs_ch2 // findRel
  publishDir params.output_dir, pattern: "*{csv,pdf}", \
             overwrite:true, mode:'copy'
  output:
      file("${base}-batch.tex")      into report_batch_report_ch
      tuple file("*.csv"), file("*pdf") into report_batch_aux_ch // need to stage
  script:
    phenotype = "pheno.phe"
    batch = "batch.phe"
    base = eigenval.baseName
    batch_col = params.batch_col
    pheno_col = params.pheno_col
    template "batchReport.py"
}
