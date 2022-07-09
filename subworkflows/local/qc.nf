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
import { FIND_SNP_EXTREME_DIFFERENTIAL_MISSINGNESS } from "../../modules/qc_processes/find_snp_extreme_differential_missingness.nf"
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

if (params.case_control) {
  ccfile = params.case_control
  Channel.fromPath(ccfile).into { cc_ch; cc2_ch }
  col    = params.case_control_col
  diffpheno = "--pheno cc.phe --pheno-name $col"
  if (params.case_control.toString().contains("s3://") || params.case_control.toString().contains("az://")) {
       println "Case control file is in the cloud so we can't check it"
  } else
  if (! file(params.case_control).exists()) {
     error("\n\nThe file <${params.case_control}> given for <params.case_control> does not exist")
    } else {
      def line
      new File(params.case_control).withReader { line = it.readLine() }
      fields = line.split()
      if (! fields.contains(params.case_control_col))
	  error("\n\nThe file <${params.case_control}> given for <params.case_control> does not have a column <${params.case_control_col}>\n")
    }

} else {
  diffpheno = ""
  col = ""
  cc_ch  = Channel.value.into("none").into { cc_ch; cc2_ch }
}


phenotype_ch = getSubChannel(params.phenotype,"pheno",params.pheno_col)
batch_ch     = getSubChannel(params.batch,"batch",params.batch_col)




/* Define the command to add for plink depending on whether sexinfo is
 * available or not.
 */

if ( nullfile.contains(params.sexinfo_available) ) {
  sexinfo = "--allow-no-sex"
  extrasexinfo = ""
  println "Sexinfo not available, command --allow-no-sex\n"
} else {
  sexinfo = ""
  extrasexinfo = "--must-have-sex"
  println "Sexinfo available command"
}


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

        checkedInput.bim_ch | GET_DUPLICATE_MARKERS

        REMOVE_DUPLICATE_SNPS(
            checkedInput.bim_ch,
            GET_DUPLICATE_MARKERS.out.duplicates_ch
        )

        //FIXME Make sure this works as expected
        def missingness = [0.01,0.03,0.05]  // this is used by one of the templates


        //TODO optional analysis of X chromossome
        if (extrasexinfo == "--must-have-sex") {

            x_analy_res_ch = REMOVE_DUPLICATE_SNPS.out.qc1_ch \
                             | GET_X.out.x_chr_ch \
                             | ANALYZE_X.out.x_analy_res_ch

        } else {

            x_analy_res_ch = Channel.fromPath("0")

        }


        REMOVE_DUPLICATE_SNPS.out.ind_miss_ch \
        | GENERATE_INDIV_MISSINGNESS_PLOT.out.report_indmisspdf_ch

        REMOVE_DUPLICATE_SNPS.out.snp_miss_ch \
        | GENERATE_SNP_MISSINGNESS_PLOT.out.report_snpmiss_ch

        qc1_ch = REMOVE_DUPLICATE_SNPS.out.qc1_ch

        qc1_ch \
        | IDENTIFY_INDIV_DISC_SEX_INFO.out.hwe_stats_ch \
        | SHOW_HWE_STATS.out.report_inithwe_ch

        qc1_ch \
        | IDENTIFY_INDIV_DISC_SEX_INFO.out.failed_sex_ch

        qc1_ch \
        | IDENTIFY_INDIV_DISC_SEX_INFO.out.batchrep_missing_ch \
        | BATCH_PROC

        //TODO Merge with the other submodules
        qc1_ch \
        | IDENTIFY_INDIV_DISC_SEX_INFO.out.failed_sex_ch \
        | REMOVE_QC_INDIVS

        qc1_ch \
        | GET_INIT_MAF.out.init_freq_ch \
        | SHOW_INIT_MAF.out.report_initmaf_ch

        qc2_ch = REMOVE_DUPLICATE_SNPS.out.qc1_ch \
                | REMOVE_QC_PHASE1.out.qc2_ch

        qc2_ch | COMP_PCA.out.pcares | DRAW_PCA

        qc2_ch | COMP_PCA.out.out_only_pcs_ch | BATCH_PROC

        qc2_ch | COMP_PCA.out.out_only_pcs_ch | PRUNE_FOR_IBD.out.find_rel_ch | BATCH_PROC


        qc2_ch \
        | COMP_PCA.out.out_only_pcs_ch \
        | PRUNE_FOR_IBD.out.find_rel_ch \
        | FIND_RELATED_INDIV.out.related_indivs_ch \
        | BATCH_PROC

        qc3_ch = REMOVE_DUPLICATE_SNPS.out.qc1_ch \
                | REMOVE_QC_PHASE1.out.qc2_ch \
                | COMP_PCA.out.out_only_pcs_ch \
                | PRUNE_FOR_IBD.out.find_rel_ch \
                | FIND_RELATED_INDIV.out.related_indivs_ch \
                | REMOVE_QC_INDIVS.out.qc3_ch

        qc3_ch \
        | CALCULATE_SNP_SKEW_STATUS.out.clean_diff_miss_plot_ch \
        | GENERATE_DIFFERENTIAL_MISSINGNESS_PLOT

        qc3_ch \
        | CALCULATE_SNP_SKEW_STATUS.out.clean_diff_miss_plot_ch \
        | FIND_SNP_EXTREME_DIFFERENTIAL_MISSINGNESS.out.skewsnps_ch

        qc3_ch \
        | CALCULATE_SNP_SKEW_STATUS.out.clean_diff_miss_plot_ch \
        | FIND_SNP_EXTREME_DIFFERENTIAL_MISSINGNESS.out.skewsnps_ch \
        | REMOVE_SKEW_SNPS.out.qc4_ch \
        | CALCULATE_MAF.out.maf_plot_ch \
        | GENERATE_MAF_PLOT

        qc3_ch \
        | CALCULATE_SNP_SKEW_STATUS.out.clean_diff_miss_plot_ch \
        | FIND_SNP_EXTREME_DIFFERENTIAL_MISSINGNESS.out.skewsnps_ch \
        | REMOVE_SKEW_SNPS.out.qc4_ch \
        | OUT_MD5


        qc2_ch \
        | CALCULATE_SAMPLE_HETEROZYGOSITY.out.hetero_check_ch \
        | GENERATE_MISS_HET_PLOT

        qc2_ch \
        | CALCULATE_SAMPLE_HETEROZYGOSITY.out.hetero_check_ch \
        | GET_BAD_INDIV_MISSING_HET.out.failed_miss_het


        qc2_ch | PRODUCE_REPORTS

        qc1_ch | PRODUCE_REPORTS


        //FIXME Make sure this is working as expected
        def mperm_header=" CHR                               SNP         EMP1         EMP2 "

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
