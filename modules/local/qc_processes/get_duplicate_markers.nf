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
        //def remove_on_bp = "??" FIX
        """
        dups.py \\
        --inpfname ${inpfname}
        --outfname ${outfname}
        --remove_on_bp ${remove_on_bp}
        """
}
