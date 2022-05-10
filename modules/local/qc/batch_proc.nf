process batchProc {
  memory plink_mem_req
  input:
    tuple path(eigenval), path(eigenvec) from pcares1
    tuple path(imiss), path(lmiss), path(sexcheck_report) from batchrep_missing_ch
    path "pheno.phe" from phenotype_ch    // staged input path
    path "batch.phe" from batch_ch        // staged input path
    path genome    from batch_rel_ch    // pruneForIBD
    path pkl       from x_analy_res_ch  // analyseX
    path rem_indivs from related_indivs_ch2 // findRel
  publishDir params.output_dir, pattern: "*{csv,pdf}", \
             overwrite:true, mode:'copy'
  output:
      path("${base}-batch.tex")      into report_batch_report_ch
      tuple path("*.csv"), path("*pdf") into report_batch_aux_ch // need to stage
  script:
    phenotype = "pheno.phe"
    batch = "batch.phe"
    base = eigenval.baseName
    batch_col = params.batch_col
    pheno_col = params.pheno_col
    template "batchReport.py"
}
