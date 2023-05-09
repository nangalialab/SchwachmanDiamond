args = commandArgs(trailingOnly=TRUE)
if(length(args)<5 || length(args)>7){
	stop("Rscript run_patient.R <patient prefix(e.g. PD6646)> <config dir> <data dir> <output dir> <MINACF> <NORMAL_ACF|optional(0.1)> <comma delimited exclude list:a,c (optional)>")
}
#source("trees_analysis.R")
source("colony_qc.R")

prefix=args[1]#"PD37580"
config_dir=args[2]
data_dir=args[3]
output_dir=sprintf("%s/%s",args[4],prefix)
if(!file.exists(output_dir)){
  dir.create(output_dir,recursive = TRUE)
}
exclude=NULL
if(file.exists(sprintf("%s/%s_qc_plots.pdf",output_dir,prefix))){

  fexc=get_exclude_file()
  if(file.exists(fexc)){
  exclude=sort(unlist(strsplit(readLines(fexc),",")))
  }
  fexc=get_tree_based_exclude_file()
  if(file.exists(fexc)){
    exclude=sort(unique(c(exclude,unlist(strsplit(readLines(fexc) %>% (function(x) x[grep("^#",x,invert = TRUE)]),",")))))
  }

  versions=as.integer(gsub("^v","",system(sprintf("ls %s | grep ^v",output_dir),intern=TRUE)))
  versions=versions[!is.na(versions)]
  if(length(versions)>0){
    lastv=max(versions)+1
  }else{
    lastv=1
  }
  cat(sprintf("archiving previous results to %s/v%d\n",output_dir,lastv))
  system(sprintf("mkdir %s/v%d;cp %s/%s_*.* %s/v%d/.",output_dir,lastv,output_dir,prefix,output_dir,lastv))
}
cat("Logging to:",sprintf("%s/%s_qc.log",output_dir,prefix),"\n")
fo=file(sprintf("%s/%s_qc.log",output_dir,prefix),open="wt")
sink(fo,split=TRUE)
sink(fo,type="message")


b_rho=FALSE
MINACF=as.numeric(args[5])
NORMAL_ACF=0.1
if(length(args)>5){
  NORMAL_ACF=as.numeric(args[6])
    if(length(args)>6){
      exclude=unlist(strsplit(args[6],","))
    }
}else{
  args=c(args,NORMAL_ACF,paste(exclude,collapse=","))
}
cat("Starting filtering & genotyping @",format(Sys.time()),"\n")
cat("Args:\n",paste(args,collapse="\n"),"\n")
inf=load_data(config_dir,data_dir,prefix,exclude)
inf$meta$NORMAL_ACF=NORMAL_ACF
inf=run_filters(inf )
inf=do_filter(inf)
inf=do_genotyping_basic_binomial(inf)
inf=do_loh(inf )
if(b_rho){
inf=apply_rho(inf,"snv")
inf=apply_rho(inf,"indel")
}
inf$meta$MINACF=MINACF
inf$meta$MINVAF=0.5*MINACF
save(inf,file=sprintf("%s/%s_qc_TMP.Rdat",output_dir,prefix))

do_qc_plots(inf,output_dir)

pdx=merge_snps_and_indels(inf)
save(inf,pdx,file=sprintf("%s/%s_qc.Rdat",output_dir,prefix))
cat("Too manually exclude write loci in the form <chr>:<pos> to",get_manual_exclude_file())
cat("Finished @",format(Sys.time()),"\n")
sink(type = "message")
sink()
