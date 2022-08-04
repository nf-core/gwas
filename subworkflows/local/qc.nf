// Import modules from modules/qc_input_validations
import { NO_SAMPLESHEET } from "../../modules/qc_input/no_samplesheet.nf"
import { SAMPLESHEET } from "../../modules/qc_input/samplesheet.nf"


// Import functions from modules/qc_utils
import { getres; checkColumnHeader; checkSampleSheet; getSubChannel } from "../../modules/qc_utils.nf"

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

workflow QC_INPUT_VALIDATION {
    main:

            if (params.idpat ==  "0")
                idpat   = "(.*)"
            else
                idpat   = params.idpat

            def K = "--keep-allele-order"

            def nullfile = [false,"False","false", "FALSE",0,"","0","null",null]


            if (nullfile.contains(params.samplesheet))
                samplesheet = "0"
            else {
                samplesheet = params.samplesheet
                checkSampleSheet(samplesheet)
            }

            def idfiles = [params.batch,params.phenotype]

            idfiles.each { checkColumnHeader(it,['FID','IID']) }

            println "The batch file is ${params.batch}"


            def max_plink_cores = params.max_plink_cores
            def plink_mem_req   = params.plink_mem_req
            def other_mem_req   = params.other_mem_req
            def pi_hat          = params.pi_hat
            def super_pi_hat    = params.super_pi_hat
            def cut_diff_miss   = params.cut_diff_miss
            def f_lo_male       = params.f_lo_male
            def f_hi_female     = params.f_hi_female
            def remove_on_bp    = params.remove_on_bp

            // def allowed_params= ["AMI","accessKey","batch","batch_col","bootStorageSize","case_control","case_control_col", "chipdescription", "cut_het_high","cut_get_low","cut_maf","cut_mind","cut_geno","cut_hwe","f_hi_female","f_lo_male","cut_diff_miss","cut_het_low", "help","input_dir","input_pat","instanceType","manifest", "maxInstances", "max_plink_cores","high_ld_regions_fname","other_mem_req","output", "output_align", "output_dir","phenotype","pheno_col","pi_hat", "plink_mem_req","region","reference","samplesheet", "scripts","secretKey","sexinfo_available", "sharedStorageMount","strandreport","work_dir","max_forks","big_time","super_pi_hat","samplesize","idpat","newpat","access-key","secret-key","instance-type","boot-storage-size","max-instances","shared-storage-mount","gemma_num_cores","remove_on_bp","queue","data","pheno","gc10"]

            // params.each { parm ->
            // if (! allowed_params.contains(parm.key)) {
            //         println "Check $parm  ************** is it a valid parameter -- are you using one rather than two - signs or vice-versa";
            //     }
            // }


            if (params.case_control) {

                ccfile = params.case_control

                Channel.fromPath(ccfile).into { cc_ch; cc2_ch }

                col    = params.case_control_col

                diffpheno = "--pheno cc.phe --pheno-name $col"

                if (params.case_control.toString().contains("s3://") || params.case_control.toString().contains("az://")) {
                    println "Case control file is in the cloud so we can't check it"
                } else if (! file(params.case_control).exists()) {
                    error("\n\nThe file <${params.case_control}> given for <params.case_control> does not exist")
                } else {
                    def line
                    new File(params.case_control).withReader { line = it.readLine() }
                    fields = line.split()
                    if (! fields.contains(params.case_control_col))
                        error("\n\nThe file <${params.case_control}> given for <params.case_control> does not have a column <${params.case_control_col}>\n")
                }
            } else {
                def diffpheno = ""
                def col = ""
                cc_ch  = Channel.value.into("none").into { cc_ch; cc2_ch }
            }


            ch_phenotype = getSubChannel(params.phenotype,"pheno",params.pheno_col)
            ch_batch     = getSubChannel(params.batch,"batch",params.batch_col)



//---- Modification of variables for pipeline -------------------------------//



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
            def inpat = "${params.input_dir}/${params.input_pat}"



            // Checks if the file exists
            checker = { fn ->
            if (fn.exists())
                return fn;
                else
                error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
            }



            if (inpat.contains("s3://") || inpat.contains("az://")) {
                print "Here"
                this_checker = { it -> return it}
                } else {
                this_checker = checker
            }


            Channel .fromFilePairs("${inpat}.{bed,bim,fam}",size:3, flat : true)
                    { file -> file.baseName }  \
                    .ifEmpty { error "No matching plink files" }        \
                    .map { a -> [this_checker(a[1]), this_checker(a[2]), this_checker(a[3])] }\
                    .multiMap {  it ->
                        raw_ch: it
                        bim_ch: it[1]
                        inpmd5ch : it
                    }.set {checked_input}


            if (samplesheet != "0")  {
                sample_sheet_ch = file(samplesheet)

                SAMPLESHEET(sample_sheet_ch)

                ch_poorgc10 = SAMPLESHEET(sample_sheet_ch).out.poorgc10_ch
                ch_report_poorgc10 = SAMPLESHEET(sample_sheet_ch).out.report_poorgc10_ch

            } else {

               ch_poorgc10 = NO_SAMPLESHEET().out.poorgc10_ch
               ch_report_poorgc10 = NO_SAMPLESHEET().out.report_poorgc10_ch

            }

    emit:
        checked_input_ch = checked_input
        poor_gc_10_ch = ch_poorgc10
        phenotype_ch = ch_phenotype
        batch_ch = ch_batch

}


workflow QC_PROCESSES {

    take:
        checked_input
        poor_gc_10_ch
        phenotype_ch
        batch_ch


    main:

        checked_input.input_md5_ch | IN_MD5.out.report_input_md5_ch

        checked_input.bim_ch | GET_DUPLICATE_MARKERS

        REMOVE_DUPLICATE_SNPS(
            checked_input.bim_ch,
            GET_DUPLICATE_MARKERS.out.duplicates_ch
        )

        qc1_ch = REMOVE_DUPLICATE_SNPS.out.qc1_ch


        //FIXME Make sure this works as expected
        def missingness = [0.01,0.03,0.05]  // this is used by one of the templates


        //TODO optional analysis of X chromosome
        if (extrasexinfo == "--must-have-sex") {

            x_analy_res_ch = qc1_ch \
                             | GET_X.out.x_chr_ch \
                             | ANALYZE_X.out.x_analy_res_ch

        } else {

            x_analy_res_ch = Channel.fromPath("0")

        }


        REMOVE_DUPLICATE_SNPS.out.ind_miss_ch \
        | GENERATE_INDIV_MISSINGNESS_PLOT.out.report_indmisspdf_ch

        REMOVE_DUPLICATE_SNPS.out.snp_miss_ch \
        | GENERATE_SNP_MISSINGNESS_PLOT.out.report_snpmiss_ch


        qc1_ch \
        | IDENTIFY_INDIV_DISC_SEX_INFO.out.hwe_stats_ch \
        | SHOW_HWE_STATS.out.report_inithwe_ch

        qc1_ch \
        | IDENTIFY_INDIV_DISC_SEX_INFO.out.failed_sex_ch

        qc1_ch \
        | IDENTIFY_INDIV_DISC_SEX_INFO.out.batchrep_missing_ch

        qc1_ch \
        | GET_INIT_MAF.out.init_freq_ch \
        | SHOW_INIT_MAF.out.report_initmaf_ch

        qc2_ch = REMOVE_DUPLICATE_SNPS.out.qc1_ch \
                | REMOVE_QC_PHASE1.out.qc2_ch

        qc2_ch | COMP_PCA.out.pcares | DRAW_PCA

        qc2_ch | COMP_PCA.out.out_only_pcs_ch | BATCH_PROC

        //TODO: The choice for find_rel_ch involves an if-else
        find_rel_ch = qc2_ch | COMP_PCA.out.out_only_pcs_ch | PRUNE_FOR_IBD.out.find_rel_ch // PRUNE_FOR_IBDLD.out.find_rel_ch


        find_rel_ch \
        | FIND_RELATED_INDIV.out.related_indivs_ch


        BATCH_PROC(
            COMP_PCA.out.pcares,
            IDENTIFY_INDIV_DISC_SEX_INFO.out.batchrep_missing_ch,
            phenotype_ch,
            batch_ch,
            find_rel_ch,
            x_analy_res_ch,
            FIND_RELATED_INDIV.out.related_indivs_ch
        )

        qc2_ch \
        | CALCULATE_SAMPLE_HETEROZYGOSITY.out.hetero_check_ch \
        | GENERATE_MISS_HET_PLOT

        qc2_ch \
        | CALCULATE_SAMPLE_HETEROZYGOSITY.out.hetero_check_ch \
        | GET_BAD_INDIV_MISSING_HET.out.failed_miss_het


        find_rel_ch | FIND_RELATED_INDIV.out.related_indivs_ch

        REMOVE_QC_INDIVS(
                        GET_BAD_INDIV_MISSING_HET.out.failed_miss_het,
                        FIND_RELATED_INDIV.out.related_indivs_ch,
                        IDENTIFY_INDIV_DISC_SEX_INFO.out.failed_sex_ch,
                        poor_gc10_ch,
                        qc2_ch
                    )


        qc3_ch = REMOVE_QC_INDIVS.out.qc3_ch

        qc3_ch \
        | CALCULATE_SNP_SKEW_STATUS.out.clean_diff_miss_plot_ch \
        | GENERATE_DIFFERENTIAL_MISSINGNESS_PLOT

        qc3_ch \
        | CALCULATE_SNP_SKEW_STATUS.out.clean_diff_miss_plot_ch \
        | FIND_SNP_EXTREME_DIFFERENTIAL_MISSINGNESS.out.skewsnps_ch


        qc4_ch = CALCULATE_SNP_SKEW_STATUS.out.clean_diff_miss_plot_ch \
                | FIND_SNP_EXTREME_DIFFERENTIAL_MISSINGNESS.out.skewsnps_ch \
                | REMOVE_SKEW_SNPS.out.qc4_ch

        qc4_ch \
        | CALCULATE_MAF.out.maf_plot_ch \
        | GENERATE_MAF_PLOT

        qc4_ch \
        | OUT_MD5


        //FIXME Make sure this is working as expected
        def mperm_header=" CHR                               SNP         EMP1         EMP2 "

    emit:
        reports_ch
}


workflow QC_WF {


    QC_INPUT_VALIDATION()

    QC_PROCESSES(
        QC_INPUT_VALIDATION.out.checked_input_ch,
        QC_INPUT_VALIDATION.out.poor_gc_10_ch,
        QC_INPUT_VALIDATION.out.phenotype_ch,
        QC_INPUT_VALIDATION.out.batch_ch,
    )

    PRODUCE_REPORTS(
        QC_PROCESSES.out.reports_ch
    )

}
