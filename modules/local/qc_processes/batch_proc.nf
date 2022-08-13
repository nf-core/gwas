process BATCH_PROC {
  publishDir params.output_dir, pattern: "*{csv,pdf}",  overwrite:true, mode:'copy'

  input:
    tuple path(eigenval), path(eigenvec)
    tuple path(imiss), path(lmiss), path(sexcheck_report)
    path("pheno.phe")    // staged input path
    path("batch.phe")        // staged input path
    path(genome)       // pruneForIBD
    path(pkl)        // analyseX
    path(rem_indivs) // findRel

  output:
    path("${base}-batch.tex"), emit: report_batch_report_ch
    tuple path("*.csv"), path("*pdf"), emit: report_batch_aux_ch // need to stage

  script:
    def phenotype = "pheno.phe"
    def batch = "batch.phe"
    def base = eigenval.baseName
    def batch_col = params.batch_col
    def pheno_col = params.pheno_col
    """
    batchReport.py \\
    --base ${base} \\
    --batch ${batch} \\
    --batch_col  ${batch_col} \\
    --phenotype ${phenotype} \\
    --pheno_col ${pheno_col} \\
    --imiss  ${imiss} \\
    --sexcheck_report ${sexcheck_report} \\
    --eigenvec  ${eigenvec} \\
    --genome  ${genome} \\
    --sx_pickle  ${pkl}
    """ 
}
