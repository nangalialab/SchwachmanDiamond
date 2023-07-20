source("load_and_annotate_tree.R")
source("project_SDS.R")
source("paper_plot_functions.R")
source("paper_export_functions.R")
source("run_timing_analysis_etc.R")
VERSION="SUBMISSION"
SELECTION_TREE_MODEL="poisson"
##  Load the data  
##  Starting point is filtered call set together with tree stored in R objects..
pdd_cache="../cache/PDD.RDS"
if(file.exists(pdd_cache)){
  cat("loading stored PDD data\n")
  PDD=readRDS(pdd_cache)
}else{
  PDD=get_PDD()  ### Relies on data with germline variants so this execution path not supported.
  cat("removing germline data...")
  PDD=lapply(PDD,function(x){x[["inf"]]=NULL;x})
  names(PDD)=sapply(PDD,function(PD) PD$INTERNAL_ID)
  cat("caching PDD to ",pdd_cache,"\n")
  saveRDS(PDD,file=pdd_cache)
}
cohort=data.frame(patient=PATIENTS,internal_id=sapply(PDD,function(PD) PD$INTERNAL_ID))
cohort$id=cohort$internal_id
cohort$N=sapply(cohort$id,function(x) length(PDD[[x]]$pdx$meta$clones))
##Alternative driver mapping names
altm=read.table("alt_driver_mappings.txt",head=TRUE,stringsAsFactors = FALSE,sep="\t")
set.seed(59666095)
##  Fit tree rate models - gets ultrametric trees and also fits per clade mutation rates (see rtreefit)
PDD=lapply(PDD,wraptreefit,niter=20000,b.fit.null = TRUE,method="nb_tree")
PDD=lapply(PDD,wraptreefit,niter=20000,b.fit.null = TRUE,method="poisson_tree")

##  report per colony mutations
mutbycolony=do.call("rbind",lapply(PDD,report_variants_per_colony))
write.table(mutbycolony,file=sprintf("../export/mutation_by_sample%s.txt",VERSION),quote=FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

##  generate per colony burden
df=do.call("rbind",lapply(PDD,function(PD) 
  PD$pdx$agedf %>% 
    filter(tip.label!="zeros")  %>% 
    mutate(internal_id=PD$INTERNAL_ID,nsub_adj=nsnv_adj.hybrid.multi) %>% 
    left_join(PD$pdx$cfg %>% mutate(tip.label=SHORT_LABEL,colony=LABEL) %>% dplyr::select(tip.label,colony) 
    )
))

df=df %>% dplyr::select(-driver) %>% left_join(altm) 
dfo=df %>% dplyr::select(colony,age_at_sample_pcy,age_at_sample_exact,internal_id,nsub_adj,meandepth,meanvaf,driver)
write.table(dfo,file=sprintf("../export/burden_by_sample_%s.txt",VERSION),quote=FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

timings=lapply(PDD[-8],function(PD) extract_timings(PD,b.use.nodes = TRUE,tree.model = "poisson_tree"))
idx=which(sapply(timings,length)>0)
timings=collate_results(timings[idx])
ds=get_driver_scheme()
timings$driver_orig=timings$driver
timings$internal_id=timings$patient
timings=timings %>% dplyr::select(-driver) %>% left_join(altm) %>% mutate(group=driver) %>% left_join(ds)
fields=c("patient","node","start","end","nchild","ncolony",	"driver_long","driver_short","driver_group","lower",	"upper",	"t_lower_lb95",	"t_lower_ub95",	"t_lower_median",	"t_upper_lb95",	"t_upper_median",	"t_upper_ub95"
         )
tmp=timings %>% mutate(driver_long=driver_orig,driver_short=driver3,driver_group=group,ncolony=N)
write.table(tmp[,fields],file="../export/selection_and_timing.txt",quote=FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

patients=c("SDS2", "SDS4" , "SDS5" , "SDS6","SDS10" )
fl=collate_driver_info(PDD[patients])
fl=fl %>% filter(upper_ub95>upper_median)
fl=fl %>% left_join(cohort[c("patient","N","internal_id")]) %>% dplyr::select(-driver) %>% left_join(altm) %>% mutate(group=driver,patient=internal_id)
fl$pct=fl$child_count/fl$N

pdf("../export/figure_3c.pdf",w=6,h=8,pointsize=10)
layout(matrix(c(1,1,2,2),ncol=2,byrow=TRUE),width=c(0.5,0.5),heights=c(0.2,0.8))
res=do_acquisition_plot_SDS_within_patient(fl %>% filter(patient=="SDS5"),b.include.stem = FALSE,b.add.bar = FALSE,b.add.axis = TRUE,xmax = 35)
res=do_acquisition_plot_SDS_merged(fl,b.include.stem = FALSE,b.add.bar = TRUE,xmax=35,i.order = -1)
dev.off()

## genes SDS and CH genes..
GENES=union(readLines("GENES_SDS.txt"),readLines("GENES_CH.txt"))
plot_all_trees_scaled(0.01,mar = c(1,2,1,3),b.fixed.scale = TRUE,d=c(21,3,12,4),figure.path=sprintf("../export/figure_2_%s.pdf",VERSION))

## export summary
lineagesummary=summarise_lineage_data(PDD,outstub = "../export/lineage_summary_CH")
## Extended figure 2..
p1=plot_comp()
p2=plot_comp(b.abs.scale = TRUE)
pdf("../export/extended_figure_2.pdf",w=8,h=6,pointsize=10)
gridExtra::grid.arrange(p1,p2,ncol=1)
dev.off()

sink("../export/sessionInfo.txt")
print(sessionInfo())
sink()