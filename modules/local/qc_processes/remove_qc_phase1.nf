process REMOVE_QC_PHASE1 {

  input:
        tuple path(bed), path(bim), path(fam)

  output:
        path("${output}*.{bed,bim,fam}"), emit: qc2_ch
        tuple path("qc1.out"), path("${output}.irem"), emit: report_qc1_ch

  script:
        def base=bed.baseName
        output = "${base}-c".replace(".","_")
        def K = "--keep-allele-order"
        if ( params.sexinfo_available != true ) {
              sexinfo = "--allow-no-sex"
              extrasexinfo = ""
            } else {
            sexinfo = ""
            extrasexinfo = "--must-have-sex"
        }

        """
        # remove really realy bad SNPs and really bad individuals
        plink $K --autosome --bfile $base $sexinfo --mind 0.1 --geno 0.1 --make-bed --out temp1
        plink $K --bfile temp1  $sexinfo --mind $params.cut_mind --make-bed --out temp2
        /bin/rm temp1.{bed,fam,bim}
        plink $K --bfile temp2  $sexinfo --geno $params.cut_geno --make-bed --out temp3
        /bin/rm temp2.{bed,fam,bim}
        plink $K --bfile temp3  $sexinfo --maf $params.cut_maf --make-bed --out temp4
        /bin/rm temp3.{bed,fam,bim}
        plink $K --bfile temp4  $sexinfo --hwe $params.cut_hwe --make-bed  --out $output
        /bin/rm temp4.{bed,fam,bim}
        cat *log > logfile
        touch tmp.irem
        cat *.irem > ${output}.irem
        qc1logextract.py logfile ${output}.irem > qc1.out
        """
}
