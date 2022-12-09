// Include modules from modules/local/qc_input_validations
include { NO_SAMPLESHEET } from "../../modules/local/qc_input_validations/no_samplesheet.nf"
include { SAMPLESHEET } from "../../modules/local/qc_input_validations/samplesheet.nf"


// Include functions from modules/local/qc_utils
include { getres } from "../../modules/local/qc_utils.nf"

// Include modules from modules/local/qc_processes
include { ANALYZE_X } from "../../modules/local/qc_processes/analyze_x.nf"
include { BATCH_PROC } from "../../modules/local/qc_processes/batch_proc.nf"
include { CALCULATE_MAF } from "../../modules/local/qc_processes/calculate_maf.nf"
include { CALCULATE_SAMPLE_HETEROZYGOSITY } from "../../modules/local/qc_processes/calculate_sample_heterozygosity.nf"
include { CALCULATE_SNP_SKEW_STATUS } from "../../modules/local/qc_processes/calculate_snp_skew_status.nf"
include { COMP_PCA } from "../../modules/local/qc_processes/comp_pca.nf"
include { CONVERT_IN_VCF } from "../../modules/local/qc_processes/convert_in_vcf.nf"
include { FIND_SNP_EXTREME_DIFFERENTIAL_MISSINGNESS } from "../../modules/local/qc_processes/find_snp_extreme_differential_missingness.nf"
include { FIND_HWE_OF_SNPS } from "../../modules/local/qc_processes/find_hwe_of_snps.nf"
include { FIND_RELATED_INDIV } from "../../modules/local/qc_processes/find_related_indiv.nf"
include { GET_BAD_INDIVS_MISSING_HET } from "../../modules/local/qc_processes/get_bad_indivs_missing_het.nf"
include { GET_DUPLICATE_MARKERS } from "../../modules/local/qc_processes/get_duplicate_markers.nf"
include { GET_INIT_MAF } from "../../modules/local/qc_processes/get_init_maf.nf"
include { GET_X } from "../../modules/local/qc_processes/get_x.nf"
include { IDENTIFY_INDIV_DISC_SEX_INFO } from "../../modules/local/qc_processes/identify_indiv_disc_sex_info.nf"
include { PRUNE_FOR_IBD } from "../../modules/local/qc_processes/prune_for_ibd.nf"
include { PRUNE_FOR_IBDLD } from "../../modules/local/qc_processes/prune_for_ibdld.nf"
include { REMOVE_DUPLICATE_SNPS } from "../../modules/local/qc_processes/remove_duplicate_snps.nf"
include { REMOVE_QC_INDIVS } from "../../modules/local/qc_processes/remove_qc_indivs.nf"
include { REMOVE_QC_PHASE1 } from "../../modules/local/qc_processes/remove_qc_phase1.nf"
include { REMOVE_SKEW_SNPS } from "../../modules/local/qc_processes/remove_skew_snps.nf"
include { SHOW_HWE_STATS } from "../../modules/local/qc_processes/show_hwe_stats.nf"
include { SHOW_INIT_MAF } from "../../modules/local/qc_processes/show_init_maf.nf"

// Include modules from modules/local/qc_report
include { DRAW_PCA } from "../../modules/local/qc_report/draw_pca.nf"
include { GENERATE_DIFFERENTIAL_MISSINGNESS_PLOT } from "../../modules/local/qc_report/generate_differential_missingness_plot.nf"
include { GENERATE_HWE_PLOT } from "../../modules/local/qc_report/generate_hwe_plot.nf"
include { GENERATE_INDIV_MISSINGNESS_PLOT } from "../../modules/local/qc_report/generate_indiv_missingness_plot.nf"
include { GENERATE_MAF_PLOT } from "../../modules/local/qc_report/generate_maf_plot.nf"
include { GENERATE_MISS_HET_PLOT } from "../../modules/local/qc_report/generate_miss_het_plot.nf"
include { GENERATE_SNP_MISSINGNESS_PLOT } from "../../modules/local/qc_report/generate_snp_missingness_plot.nf"
include { IN_MD5 } from "../../modules/local/qc_report/in_md5.nf"
include { OUT_MD5 } from "../../modules/local/qc_report/out_md5.nf"
include { PRODUCE_REPORTS } from "../../modules/local/qc_report/produce_reports.nf"



//======================
// Definitions of some utils
//======================

// This method first checks that the data file has the stated column
// If so, it creates a channel for it
// NB: if the file is in S3 we cannot do the test since Groovy does not
// allow us to access the file directly
def getSubChannel = { parm, parm_name, col_name ->
    if (parm.toString().contains("s3://")) {
        println "The file <$parm> is in S3 so we cannot do a pre-check";
        return Channel.fromPath(parm);
    }
    if (parm.toString().contains("az://")) {
        println "The file <$parm> is in Azure so we cannot do a pre-check";
        return Channel.fromPath(parm);
    }
    if ((parm==0) || (parm=="0") || (parm==false) || (parm=="false")) {
        filename = "emptyZ0${parm_name}.txt";
        new File(filename).createNewFile()
        new_ch = Channel.fromPath(filename);

    } else {
        if (! file(parm).exists()) {
        error("\n\nThe file <$parm> given for <params.${parm_name}> does not exist")
        } else {
        def line
        new File(parm).withReader { line = it.readLine() }
        fields = line.split()
        if (! fields.contains(col_name))
        error("\n\nThe file <$parm> given for <params.${parm_name}> does not have a column <${col_name}>\n")
        }
        new_ch = Channel.fromPath(parm);
    }
    return new_ch;
}


def getConfig = {
    all_files = workflow.configFiles.unique()
    text = ""
    all_files.each { fname ->
        base = fname.baseName
        curr = "\n\n*-subsection{*-protect*-url{$base}}@.@@.@*-footnotesize@.@*-begin{verbatim}"
        file(fname).eachLine { String line ->
        if (line.contains("secretKey")) { line = "secretKey='*******'" }
            if (line.contains("accessKey")) { line = "accessKey='*******'" }
            curr = curr + "@.@"+line
        }
        curr = curr +"@.@*-end{verbatim}\n"
        text = text+curr
    }
    return text
}

def checkColumnHeader(fname, columns, nullfile) {
    if (workflow.profile == "awsbatch") return;
    if (fname.toString().contains("s3://")) return;
    //FIXME The presence of nullfile (undeclared function param)
    if (nullfile.contains(fname)) return;
    new File(fname).withReader { line = it.readLine().tokenize() }
    problem = false;
    columns.each { col ->
        if (! line.contains(col) ) {
        println "The file <$fname> does not contain the column <$col>";
        problem=true;
        }
        if (problem)
        System.exit(2)
    }
}


def checkSampleSheet(fname)  {
    if (workflow.profile == "awsbatch") return;
    if (fname.contains("s3://") )return;
    if (fname.contains("az://") ) return;
    //FIXME The presence of nullfile (undeclared function param)
    if (nullfile.contains(fname) || fname.contains(".xls")) return;
    new File(fname).withReader { line = it.readLine()}
    problem  = false
    prob_str = ""
    if (! line.contains(",")) {
        problem = true;
        prob_str = "If given as a CSV file, it must be comma-separated\n";
    }
    headers = line.tokenize(",")
    headers.each { println it}
    if (!(headers.contains("Institute Sample Label") ||
        (headers.contains("Sample Plate") && headers.contains("Well")))) {
        problem= true
        prob_str = prob_str + "Column headers must include 'Institute Sample Label'  or both 'Sample Plate' and 'Well'"
    }
    if (problem)  {
        println "There's a problem with the sample sheet <$fname>."
        println prob_str;
        //FIXME
        System.exit(1)
    }
}


//======================
// Definitions of local workflows
//======================

workflow QC_INPUT_VALIDATION {
    main:

            //======================
            // Validation checks and channel construction
            //======================


            if (params.idpat ==  "0")
                idpat   = "(.*)"
            else
                idpat   = params.idpat

            def K = "--keep-allele-order"


            //----------------------------
            // Adapt as per the samplesheet
            //----------------------------

            def nullfile = [false,"False","false", "FALSE",0,"","0","null",null]

            if (nullfile.contains(params.samplesheet))
                samplesheet = "0"
            else {
                samplesheet = params.samplesheet
                checkSampleSheet(samplesheet)
            }



            //----------------------------

            def idfiles = [params.batch,params.phenotype]

            //FIXME Disabled the check to avoid the compilation/scope issues
            // idfiles.each { checkColumnHeader(it, ['FID','IID', nullfile]) }

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

                cc_ch = Channel.fromPath(ccfile)

                def col    = params.case_control_col

                def diffpheno = "--pheno cc.phe --pheno-name $col"

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
                cc_ch  = Channel.value("none")
            }


            ch_phenotype = getSubChannel(params.phenotype,"pheno",params.pheno_col)
            ch_batch     = getSubChannel(params.batch,"batch",params.batch_col)



//---- Modification of variables for pipeline -------------------------------//



            // Define the command to add for plink depending on whether sexinfo is available or not.

            if ( nullfile.contains(params.sexinfo_available) ) {
                sexinfo = "--allow-no-sex"
                extrasexinfo = ""
                println "Sexinfo not available, command --allow-no-sex\n"
            } else {
                sexinfo = ""
                extrasexinfo = "--must-have-sex"
                println "Sexinfo available command"
            }


            // Get the input files -- could be a glob
            // We match the bed, bim, fam file -- order determined lexicographically
            // not by order given, we check that they exist and then
            // send the all the files to raw_ch and just the bim file to bim_ch
            def inpat = "${params.input_dir}/${params.input_pat}"




            // Checks if the file exists
            if (inpat.contains("s3://") || inpat.contains("az://")) {
                    def this_checker = { it -> return it}
                } else {
                    def this_checker = { fn ->
                                        if (fn.exists())
                                            return fn;
                                            else
                                            error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
                                        }
            }


            Channel .fromFilePairs("${inpat}.{bed,bim,fam}",size:3, flat : true)
                    { file -> file.baseName }  \
                    .ifEmpty { error "No matching plink files" }        \
                    // .map { a -> [this_checker(a[1]), this_checker(a[2]), this_checker(a[3])] } \
                    .multiMap {  it ->
                        raw_ch: it
                        bim_ch: it[2]
                        inpmd5ch : it
                    }.set {checked_input}


            if (samplesheet != "0")  {
                sample_sheet_ch = file(samplesheet)

                SAMPLESHEET(sample_sheet_ch)

                ch_poorgc10 = SAMPLESHEET.out.poorgc10_ch
                ch_report_poorgc10 = SAMPLESHEET.out.report_poorgc10_ch

            } else {

                NO_SAMPLESHEET()

                ch_poorgc10 = NO_SAMPLESHEET.out.poorgc10_ch
                ch_report_poorgc10 = NO_SAMPLESHEET.out.report_poorgc10_ch

            }

    emit:
        checked_input_md5_ch = checked_input.inpmd5ch
        checked_input_bim_ch = checked_input.bim_ch
        poor_gc_10_ch = ch_poorgc10
        phenotype_ch = ch_phenotype
        batch_ch = ch_batch

}


workflow QC_PROCESSES {

    take:
        checked_input_md5_ch
        checked_input_bim_ch
        poor_gc_10_ch
        phenotype_ch
        batch_ch


    main:

        //checked_input_md5_ch | IN_MD5.out.report_input_md5_ch
        IN_MD5(checked_input_md5_ch)
        
        checked_input_bim_ch | GET_DUPLICATE_MARKERS
        
        REMOVE_DUPLICATE_SNPS(
            checked_input_md5_ch,
            GET_DUPLICATE_MARKERS.out.duplicates_ch
        )
        
        qc1_ch = REMOVE_DUPLICATE_SNPS.out.qc1_ch

        /*
        //FIXME Make sure this works as expected
        def missingness = [0.01,0.03,0.05]  // this is used by one of the templates


        //TODO optional analysis of X chromosome
        if (extrasexinfo == "--must-have-sex") {

            x_analy_res_ch = qc1_ch \
                             | GET_X.out.x_chr_ch \
                             | ANALYZE_X.out.x_analy_res_ch

        } else {

            x_analy_res_ch = Channel.fromPath("0")

        
        */

        REMOVE_DUPLICATE_SNPS.out.ind_miss_ch | GENERATE_INDIV_MISSINGNESS_PLOT
        //.out.report_indmisspdf_ch
        
        REMOVE_DUPLICATE_SNPS.out.snp_miss_ch \
        | GENERATE_SNP_MISSINGNESS_PLOT
        //.out.report_snpmiss_ch

        
        qc1_ch | IDENTIFY_INDIV_DISC_SEX_INFO
        
        IDENTIFY_INDIV_DISC_SEX_INFO.out.hwe_stats_ch 
        | SHOW_HWE_STATS //.out.report_inithwe_ch
    
        
        //.out.failed_sex_ch

        qc1_ch | GET_INIT_MAF
        
        GET_INIT_MAF.out.init_freq_ch | SHOW_INIT_MAF
        
        REMOVE_DUPLICATE_SNPS.out.qc1_ch | REMOVE_QC_PHASE1
        
        qc2_ch = REMOVE_QC_PHASE1.out.qc2_ch

        qc2_ch | COMP_PCA

        // DRAW_PCA(COMP_PCA.out.pcares, case_control)
        
        // BATCH_PROC() <- 7 input channels
        /*
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
        | GET_BAD_INDIVS_MISSING_HET.out.failed_miss_het


        find_rel_ch | FIND_RELATED_INDIV.out.related_indivs_ch

        REMOVE_QC_INDIVS(
                        GET_BAD_INDIVS_MISSING_HET.out.failed_miss_het,
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
    */    
}


workflow QC_WF {


    QC_INPUT_VALIDATION()

    QC_PROCESSES(
         QC_INPUT_VALIDATION.out.checked_input_md5_ch,
         QC_INPUT_VALIDATION.out.checked_input_bim_ch,
         QC_INPUT_VALIDATION.out.poor_gc_10_ch,
         QC_INPUT_VALIDATION.out.phenotype_ch,
         QC_INPUT_VALIDATION.out.batch_ch
    )

    // PRODUCE_REPORTS(
    //     QC_PROCESSES.out.reports_ch
    //
} 
