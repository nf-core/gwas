import REMOVE_DUPLICATE_SNPS ...
import GET_X ...

workflow QC_PROCESSES {

    REMOVE_DUPLICATE_SNPS(input_TODO)

    GET_X(
        REMOVE_DUPLICATE_SNPS.out.qc1_ch
    )


    ANALYZE_X(
        GET_X.out.x_chr_ch
    )


}

workflow QC_INPUT_VALIDATION {

}


workflow QC_WF {

    QC_INPUT_VALIDATION()

    QC_PROCESSES()

    PRODUCE_REPORTS(
        QC_PROCESSES.out.reports_ch
    )

}
