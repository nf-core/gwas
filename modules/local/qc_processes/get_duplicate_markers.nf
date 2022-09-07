/* Process to find duplicates. *
 * Inputs:
 * - bim: the bim file
 * Outputs:
 * - duplicates.snps    : A possibly empty file with a list of SNPs
 */
process GET_DUPLICATE_MARKERS {

  input:
        path(inpfname)

  output:
        path("${base}.dups"), emit: duplicates_ch

  script:
        def base     = inpfname.baseName
        def outfname = "${base}.dups"

        template "dups.py"
}