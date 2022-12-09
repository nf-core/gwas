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
        path("${prefix}.dups"), emit: duplicates_ch

  script:

        prefix     = inpfname.baseName
        def outfname = "${prefix}.dups"
        //def remove_on_bp = "??" FIX
        """
        dups.py \\
        --inpfname ${inpfname} \\
        --outfname ${outfname} \\
        --remove_on_bp ${params.remove_on_bp}
        """
}
