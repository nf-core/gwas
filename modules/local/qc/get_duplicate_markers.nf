/* Process to find duplicates. *
 * Inputs:
 * - bim: the bim file
 * Outputs:
 * - duplicates.snps    : A possibly empty file with a list of SNPs
 */
process getDuplicateMarkers {
  memory other_mem_req
  publishDir params.output_dir, pattern: "*dups", \
             overwrite:true, mode:'copy'
  input:
    path(inpfname) from checked_input.bim_ch
  output:
    path("${base}.dups") into duplicates_ch
  script:
     base     = inpfname.baseName
     outfname = "${base}.dups"
     template "dups.py"
}
