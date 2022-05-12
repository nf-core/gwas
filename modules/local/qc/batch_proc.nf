process batchProc {
  memory plink_mem_req

  input:
    tuple path(eigenval), path(eigenvec)
    tuple path(imiss), path(lmiss), path(sexcheck_report)
    path("pheno.phe")    // staged input path
    path("batch.phe")        // staged input path
    path(genome)       // pruneForIBD
    path(pkl)        // analyseX
    path(rem_indivs) // findRel

  publishDir params.output_dir, pattern: "*{csv,pdf}", \
             overwrite:true, mode:'copy'
  output:
    path("${base}-batch.tex"),emit: report_batch_report_ch
    tuple path("*.csv"), path("*pdf"), emit: report_batch_aux_ch // need to stage

  script:
    phenotype = "pheno.phe"
    batch = "batch.phe"
    base = eigenval.baseName
    batch_col = params.batch_col
    pheno_col = params.pheno_col
    template "batchReport.py"
}
