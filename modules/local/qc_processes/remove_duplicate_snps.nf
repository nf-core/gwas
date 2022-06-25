/*  Process to remove duplicate SNPs.
 * Inputs
 *   -  raw files from from user-specified data
 *   -  list of duplicates comes from getDuplicateMarkers
 * Outputs:
 *   nodups.{bed,bim,fam} (PLINK file without duplicates) and
 *   qc.log log file
 */

process REMOVE_DUPLICATE_SNPS {
  memory plink_mem_req

  input:
   tuple path(bed), path(bim), path(fam)
   path(dups)

  output:
    tuple  path("${nodup}.bed"), path("${nodup}.bim"), path("${nodup}.fam"), emit: qc1_ch
    tuple file("${base}.orig"), path(dups), emit: report_dups_ch
    path("${nodup}.lmiss"), emit: snp_miss_ch
    path("${nodup}.imiss"), emit: ind_miss_ch1

  script:
   def base    = bed.baseName
   def nodup   = "${base}-nd"
   """
    plink $K --bfile $base $sexinfo $extrasexinfo --exclude $dups --missing --make-bed --out $nodup
    wc -l ${base}.bim > ${base}.orig
    wc -l ${base}.fam >> ${base}.orig
   """
}
