/*  Process to remove duplicate SNPs.
 * Inputs
 *   -  raw files from from user-specified data
 *   -  list of duplicates comes from getDuplicateMarkers
 * Outputs:
 *   nodups.{bed,bim,fam} (PLINK file without duplicates) and
 *   qc.log log file
 */

process removeDuplicateSNPs {
  memory plink_mem_req
  input:
   tuple file(bed), file(bim), file(fam) from checked_input.raw_ch
   file(dups) from  duplicates_ch

  output:
    tuple  file("${nodup}.bed"),file("${nodup}.bim"),file("${nodup}.fam")\
    into (qc1_ch,qc1B_ch,qc1C_ch,qc1D_ch,qc1E_ch)
    tuple file("${base}.orig"), file(dups) into report_dups_ch
    file ("${nodup}.lmiss") into snp_miss_ch
    file ("${nodup}.imiss") into (ind_miss_ch1, ind_miss_ch2)
  script:
   base    = bed.baseName
   nodup   = "${base}-nd"
   """
    plink $K --bfile $base $sexinfo $extrasexinfo --exclude $dups --missing --make-bed --out $nodup
    wc -l ${base}.bim > ${base}.orig
    wc -l ${base}.fam >> ${base}.orig
   """
}
