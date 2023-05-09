PATIENTS=c("PD43294",  "PD43975",  "PD43295",  "PD43293",  "PDHW1563" ,"PD45887", "PD44706" , "PD42191"  ,"PD45889", "PD45886" )
CV=c("missense","nonsense","ess_splice","frameshift","inframe","loh","start_lost","cna","stop_lost")
GENES=readLines("GENES.txt")
get_driver_scheme2=function(){
  driver.scheme=read.table("driver_scheme_simple.txt",head=T,stringsAsFactors = FALSE,sep="\t")
  n=max(driver.scheme$number)
  pallete=c(RColorBrewer::brewer.pal(9,"Set1")[-6],RColorBrewer::brewer.pal(8,"Dark2"))
  driver.scheme$colour=pallete[driver.scheme$number]
  driver.scheme
}

##Import age of sample data
##Add all relevant meta
add_agedf_SDS=function(PD){
  dat=read.table("../data/sample_info.txt",head=TRUE,stringsAsFactors = FALSE,sep="\t")
  tree=PD$pdx$tree_ml
  sample=PD$pdx$cfg$LABEL[match(PD$pdx$tree_ml$tip.label,PD$pdx$cfg$SHORT_LABEL)]
  if(length(which(is.na(sample)==1))){
    sample[is.na(sample)]="zeros"
  }else{
    stop("unexpected lookup error!")
  }
  ##Find fuzzy match to PDID
  matches=lapply(dat$PDID,function(x) grep(x,sample))
  idx=which(sapply(matches,length)>0)
  agedf=data.frame(tip.label=tree$tip.label,sample=sample,age_at_sample_exact=rep(1e-6,length(tree$tip.label)))
  for(i in idx){
    cat("mapping",dat$PDID[i],"to tree..")
    agedf$age_at_sample_exact[matches[[i]]]=dat$Age_years[i]
    PD$INTERNAL_ID=dat$SDSID[i] ##dat$InternalID[i]
    PD$LONG_DESC=sprintf("%s: %s",dat$SDSID[i],dat$Germline_SBDS_genotype[i])
  }
  agedf$age_at_sample_pcy=agedf$age_at_sample_exact+(MEAN_AGE_AT_DELIVERY_DAYS)/365.25
  agedf$age_at_sample_pcy[agedf$tip.label=="zeros"]=1e-6
  PD$pdx$agedf=agedf
  PD
}

add_agedf=add_agedf_SDS

get_PDD=function(){
  PDD=lapply(PATIENTS,
             function(patient){
               if(patient=="PD42191"){
                 PD=get_pd(patient,b.remove.germline = FALSE)
                 ## Remove copy number events as the list is too long but also not exhaustive for this patient.
                 PD$pdx$meta$LOH=list()
                 PD$pdx$meta$CNA=list()
                 PD$nodes=PD$nodes %>% dplyr::select(node,driver,status,driver2,driver3,child_count)
                 PD
               }else{
                 PD=get_pd(patient,b.remove.germline = TRUE)
                 PD$nodes=PD$nodes %>% dplyr::select(node,driver,status,driver2,driver3,child_count)
                 PD
               }}
  )
  PDD=lapply(PDD,add_adjustment_models)
  ##PDD=lapply(PATIENTS,get_PD)
  names(PDD)=PATIENTS
  PDD
}

add_derived_fields2=function(PD){
  nodes=PD$nodes %>% dplyr::select(node,driver,status,driver2,driver3,child_count) %>% mutate(status=1)
  tree=get_colored_markup_tree2(PD$pdx$tree_ml,nodes)
  nodes=rbind(data.frame(node=-1,driver="WT",status=1,driver2="WT",driver3="WT",child_count=-1,stringsAsFactors = FALSE),nodes)
  M=dim(nodes)[1]
  tipnode=tree$rate[match(1:dim(PD$pdx$agedf)[1],tree$edge[,2])]
  tipdesc=nodes$driver[tipnode]
  PD$pdx$agedf$driver=tipdesc
  PD$pdx$agedf$driver3=nodes$driver3[tipnode]
  PD
}


collate_driver_info=function(PDD,treemodel="poisson_tree",b.is.null=TRUE){
  inf1=do.call("rbind",lapply(PDD,function(x){
    altnull=ifelse(b.is.null,"nullmodel","altmodel")
    out=x$fit[[treemodel]][[altnull]]$summary %>% mutate(patient=x$patient)
    #out=out %>% left_join(add_cumulative_mutcount(x,x$fit[[treemodel]][[type]]$summary$node))
    idx=match("patient",colnames(out))
    cbind(data.frame(patient=out[,idx],stringsAsFactors = FALSE),out[,-idx])
    ##Add
  }
  ))
  ##Make driver correspond to events on the current node
  inf1$driver_full=inf1$driver
  #inf1$driver=gsub(":.*","",inf1$driver)
  inf1$driver=inf1$driver3
  #The driver scheme maps specific mutations to drivers..
  ds=get_driver_scheme()
  inf1$col=ds$colour[match(inf1$driver,ds$group)]
  pti=do.call("rbind",lapply(PDD,function(PD) {
    tmp=unique(PD$pdx$agedf$age_at_sample_pcy[PD$pdx$agedf$age_at_sample_pcy>0.1])
    zz=data.frame(age_at_sample_pcy=tmp)
    zz$patient=PD$patient
    zz
  }))
  pti=pti %>% group_by(patient) %>% summarise(max_age_at_sample=max(age_at_sample_pcy))
  out=inf1 %>% left_join(pti)
  out
}


shade_between=function(x,y1,y2,color){
 # browser()
  polygon(c(x, rev(x), x[1]), c(y1, rev(y2), y1[1]),
          col = color,border = NA)
}

add_within_patient_layout=function(tmp,b.cluster.top=FALSE){
  patients=unique(tmp$patient)
  tmp$pct=tmp$child_count/tmp$N
  ##  Give binomial proprotions 
  ci=t(sapply(1:length(tmp$pct),function(i) binom.test(tmp$child_count[i],tmp$N[i])$conf.int))
  tmp$pct_ub95=ci[,2]
  tmp$pct_lb95=ci[,1]
  
  tmp=tmp[order(tmp$max_age_at_sample,tmp$patient,tmp$upper_median),]
  tmp=tmp %>% mutate(cpct=NA,e1=NA,eb=NA,et=NA,d1=NA,db=NA,dt=NA,eb1=NA,et1=NA)
  for(patient in patients){
    idx=which(tmp$patient==patient)
    #browser()
    tmp$cpct[idx]=cumsum(tmp$pct[idx])
    ## different layout options here e,d,f
    # Topn justified stacked.
    tmp$e1[idx]=1-(tmp$cpct[idx]-0.5*tmp$pct[idx])
    tmp$eb1[idx]=1-tmp$cpct[idx]
    tmp$et1[idx]=1-(tmp$cpct[idx]-tmp$pct[idx])
    tmp$e1[idx]=1-((1:length(idx)))/(length(idx)+1)
    if(!b.cluster.top){
    ##  Also maximally spread out..  
    pad=(1-sum(tmp$pct[idx]))/(length(idx)+1)
    tmp$e1[idx]=1-(tmp$cpct[idx]+pad*(1:length(idx))-0.5*tmp$pct[idx])
    tmp$eb[idx]=1-(tmp$cpct[idx]+pad*(1:length(idx)))
    tmp$et[idx]=1-(tmp$cpct[idx]+pad*(1:length(idx))-tmp$pct[idx])
    }else{
      tmp$eb=tmp$eb1
      tmp$et=tmp$et1
    }
  }
  tmp
}

get_SDS_pallete=function(){
  mycol = brewer.pal(n=9,name="Set1")
  extra=brewer.pal(n=8,name="Dark2")[6:7]
  paired=brewer.pal(n=11,name="Paired")
  extra=c(paired[5],brewer.pal(n=8,name="Dark2")[8],RColorBrewer::brewer.pal(12,"Set3")[11])
  ## variation on colour scheme
  pallete = c(mycol[c(2)],"#a6cee3","#33a02c","#b2df8a", mycol[c(5)],"#8761B7", "#cab2d6",mycol[c(1,7)], "light grey",extra)
  pallete
}

get_driver_scheme=function(){
  driver.scheme=read.table("driver_scheme_simple_SDS.txt",head=T,stringsAsFactors = FALSE,sep="\t")
  pallete=get_SDS_pallete()
  n=max(driver.scheme$number)
  #pallete=c(RColorBrewer::brewer.pal(9,"Set1")[-6],RColorBrewer::brewer.pal(8,"Dark2"))
  driver.scheme$colour=toupper(pallete[driver.scheme$number])
  ch=readLines("GENES_CH.txt")
  ch=ch[which(!(ch %in% driver.scheme$group))]
  driver.scheme=rbind(driver.scheme,data.frame(group=ch,number=14,colour="black"))
}
get_driver_scheme2=get_driver_scheme

