#!/usr/bin/env nextflow

import java.nio.file.Paths;
import sun.nio.fs.UnixPath;
import java.security.MessageDigest;


filescript=file(workflow.scriptFile)
projectdir="${filescript.getParent()}"
dummy_dir="${projectdir}/../../qc/input"


// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}

def helps = [ 'help' : 'help' ]
allowed_params = ['file_gwas', 'file_ref_gzip', "output_dir","output", "input_dir", "input_pat"]
allowed_header = ['head_pval', 'head_freq', 'head_bp', 'head_chr', 'head_rs', 'head_beta', 'head_se', 'head_A1', 'head_A2', 'sep']
allowed_headnewernew = ['headnew_pval', 'headnew_freq', 'headnew_bp', 'headnew_chr', 'headnew_rs', 'headnew_beta', 'headnew_se', 'headnew_A1', 'headnew_A2']
//chro_ps 0 --bp_ps 1 --rs_ps 
allowed_posfilref=['poshead_chro_inforef', 'poshead_bp_inforef','poshead_rs_inforef']
allowed_params+=allowed_header
allowed_params+=allowed_headnewernew

if(params.file_gwas==""){
error('params.file_gwas: file contains gwas not found')
}

if(params.headnew_pval=="")params.headnew_pval=params.head_pval 
if(params.headnew_freq=="")params.headnew_freq=params.head_freq 
if(params.headnew_bp=="")params.headnew_bp=params.head_bp 
if(params.headnew_chr=="")params.headnew_chr=params.head_chr 
if(params.headnew_beta=="")params.headnew_beta=params.head_beta 
if(params.headnew_se=="")params.headnew_se=params.head_se 
if(params.headnew_A1=="")params.headnew_A1=params.head_A1 
if(params.headnew_A2=="")params.headnew_A2=params.head_A2 
//if(params.headnew_N=="")params.headnew_N=params.head_N 


gwas_chrolist = Channel.fromPath(params.file_gwas)
gwas_chrolist_ext = Channel.fromPath(params.file_gwas)
gwas_chrolist_ext = Channel.fromPath(params.file_gwas)


getListeChro(
    gwas_chrolist
).out.chrolist

chrolist2=Channel.create()
chrolist.flatMap { list_str -> list_str.readLines()[0].split() }.set { chrolist2 }

ExtractChroGWAS(
        file(gwas) from gwas_chrolist_ext
        each chro from  chrolist2
}.out.gwas_format_chro

if(params.file_ref_gzip==""){
    error('params.file_ref_gzip : file contains information for rs notnot found')
}

gwas_format_chro_rs=gwas_format_chro.combine(Channel.fromPath(params.file_ref_gzip))

ExtractRsIDChro(
    gwas_format_chro_rs
).out.rsinfo_chro

if(params.input_dir!="" || params.input_pat!=''){
print("used plink file")
bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")


}else{
    bed=file('${dummy_dir}/00')
    bim=file('${dummy_dir}/01')
    fam=file('${dummy_dir}/02')
}

fileplk= Channel.create()
Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [a[0], a[1], a[2]] }
    .set { fileplk }

rsinfo_chroplk=rsinfo_chro.combine(fileplk)

process MergeRsGwasChro{
    rsinfo_chroplk
).out.gwas_rsmerge

gwas_rsmerge_all=gwas_rsmerge.collect()

MergeAll(
    gwas_rsmerge_all
}.out



