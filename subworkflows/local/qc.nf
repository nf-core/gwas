// Import modules from modules/qc_input_validations
import { NO_SAMPLESHEET } from "../../modules/qc_input/no_samplesheet.nf"
import { SAMPLESHEET } from "../../modules/qc_input/samplesheet.nf"

// Import modules from modules/qc_processes
import { ANALYZE_X } from "../../modules/qc_processes/analyze_x.nf"
import { BATCH_PROC } from "../../modules/qc_processes/batch_proc.nf"
import { CALCULATE_MAF } from "../../modules/qc_processes/calculate_maf.nf"
import { CALCULATE_SAMPLE_HETEROZYGOSITY } from "../../modules/qc_processes/calculate_sample_heterozygosity.nf"
import { CALCULATE_SNP_SKEW_STATUS } from "../../modules/qc_processes/calculate_snp_skew_status.nf"
import { COMP_PCA } from "../../modules/qc_processes/comp_pca.nf"
import { CONVERT_IN_VCF } from "../../modules/qc_processes/convert_in_vcf.nf"
import { FIND_EXTREME_DIFFERENTIAL_MISSINGNESS } from "../../modules/qc_processes/find_extreme_differential_missingness.nf"
import { FIND_HWE_OF_SNPS } from "../../modules/qc_processes/find_hwe_of_snps.nf"
import { FIND_RELATED_INDIV } from "../../modules/qc_processes/find_related_indiv.nf"
import { GET_BAD_INDIV_MISSING_HET } from "../../modules/qc_processes/get_bad_indiv_missing_het.nf"
import { GET_DUPLICATE_MARKERS } from "../../modules/qc_processes/get_duplicate_markers.nf"
import { GET_INIT_MAF } from "../../modules/qc_processes/get_init_maf.nf"
import { GET_X } from "../../modules/qc_processes/get_x.nf"
import { IDENTIFY_INDIV_DISC_SEX_INFO } from "../../modules/qc_processes/identify_indiv_disc_sex_info.nf"
import { PRUNE_FOR_IBD } from "../../modules/qc_processes/prune_for_ibd.nf"
import { PRUNE_FOR_IBDLD } from "../../modules/qc_processes/prune_for_ibdld.nf"
import { REMOVE_DUPLICATE_SNPS } from "../../modules/qc_processes/remove_duplicate_snps.nf"
import { REMOVE_QC_INDIVS } from "../../modules/qc_processes/remove_qc_indivs.nf"
import { REMOVE_QC_PHASE1 } from "../../modules/qc_processes/remove_qc_phase1.nf"
import { REMOVE_SKEW_SNPS } from "../../modules/qc_processes/remove_skew_snps.nf"
import { SAMPLE_SHEET } from "../../modules/qc_processes/sample_sheet.nf"
import { SHOW_HWE_STATS } from "../../modules/qc_processes/show_hwe_stats.nf"
import { SHOW_INIT_MAF } from "../../modules/qc_processes/show_init_maf.nf"

// Import modules from modules/qc_report

import { DRAW_PCA } from "../../modules/qc_report/draw_pca.nf"
import { GENERATE_DIFFERENTIAL_MISSINGNESS_PLOT } from "../../modules/qc_report/generate_differential_missingness_plot.nf"
import { GENERATE_HWE_PLOT } from "../../modules/qc_report/generate_hwe_plot.nf"
import { GENERATE_INDIV_MISSINGNESS_PLOT } from "../../modules/qc_report/generate_indiv_missingness_plot.nf"
import { GENERATE_MAF_PLOT } from "../../modules/qc_report/generate_maf_plot.nf"
import { GENERATE_MISS_HET_PLOT } from "../../modules/qc_report/generate_miss_het_plot.nf"
import { GENERATE_SNP_MISSINGNESS_PLOT } from "../../modules/qc_report/generate_snp_missingness_plot.nf"
import { IN_MD5 } from "../../modules/qc_report/in_md5.nf"
import { OUT_MD5 } from "../../modules/qc_report/out_md5.nf"
import { PRODUCE_REPORTS } from "../../modules/qc_report/produce_reports.nf"


//======================
// Definitions of local workflows
//======================

//FIXME
workflow QC_INPUT_VALIDATION {

/* Get the input files -- could be a glob
 * We match the bed, bim, fam file -- order determined lexicographically
 * not by order given, we check that they exist and then
 * send the all the files to raw_ch and just the bim file to bim_ch */
inpat = "${params.input_dir}/${params.input_pat}"



if (inpat.contains("s3://") || inpat.contains("az://")) {
  print "Here"
  this_checker = { it -> return it}
} else {
  this_checker = checker
}


    Channel
   .fromFilePairs("${inpat}.{bed,bim,fam}",size:3, flat : true)
    { file -> file.baseName }  \
      .ifEmpty { error "No matching plink files" }        \
      .map { a -> [this_checker(a[1]), this_checker(a[2]), this_checker(a[3])] }\
      .multiMap {  it ->
         raw_ch: it
         bim_ch: it[1]
         inpmd5ch : it
       }.set {checked_input}

}


workflow QC_PROCESSES {

    take:
        checkedInput


    main:
        GET_DUPLICATE_MARKERS( checkedInput.bim_ch )

        REMOVE_DUPLICATE_SNPS(
            checkedInput.bim_ch,
            GET_DUPLICATE_MARKERS.out.duplicates_ch
        )

        //TODO optional analysis of X chromossome
        if (extrasexinfo == "--must-have-sex") {
            GET_X(REMOVE_DUPLICATE_SNPS.out.qc1_ch)
            ANALYZE_X(GET_X.OUT.X_chr_ch)
        )

        IDENTIFY_INDIV_DISC_SEXINFO(
            REMOVE_DUPLICATE_SNPS.out.qc1_ch
        )

        GET_INIT_MAF(
            REMOVE_DUPLICATE_SNPS.out.qc1_ch
        )


        //PLOTS
        GENERATE_SNP_MISSINGNESS_PLOT(
            REMOVE_DUPLICATE_SNPS.out.snp_miss_ch
        )

        GENERATE_INDIV_MISSINGNESS_PLOT(
            REMOVE_DUPLICATE_SNPS.out.ind_miss_ch
        )

        SHOW_INIT_MAF(
            GET_INIT_MAF.out.init_freq_ch
        )

        SHOW_HWE_STATS(
            IDENTIFY_INDIV_DISC_SEXINFO.out.hwe_stats_ch
        )


}


workflow QC_WF {


    QC_INPUT_VALIDATION()

    QC_PROCESSES(
        QC_INPUT_VALIDATION.out.checkedInput
    )

    PRODUCE_REPORTS(
        QC_PROCESSES.out.reports_ch
    )

}
