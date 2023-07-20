
get_lineages=function(tree,threshold,min.count=2){
  nh=nodeHeights(tree)
  idx=which(nh[,1]<threshold & nh[,2]>=threshold)
  nodes=tree$edge[idx,2]
  out=data.frame(node=nodes,start=nh[idx,1],end=nh[idx,2])
  if(length(nodes)>0){
    out$N=sapply(nodes,function(node) length(get_samples_in_clade(node,tree)))
  }else{
    out$N=integer()
  }
  out %>% filter(N>1)
}

get_lineages2=function(tree,nodes){
  nh=nodeHeights(tree)
  idx=match(nodes,tree$edge[,2])
  if(any(is.na(idx))){
    browser()
    stop("bad node list")
  }
  nodes=tree$edge[idx,2]
  out=data.frame(node=nodes,start=nh[idx,1],end=nh[idx,2])
  if(length(nodes)>0){
    out$N=sapply(nodes,function(node) length(get_samples_in_clade(node,tree)))
  }else{
    out$N=integer()
  }
  out %>% filter(N>1)
}
get_lineages3=function(tree,threshold,min.count=2){
  nh=nodeHeights(tree)
  idx=which(nh[,1]<threshold & nh[,2]>=threshold)
  nodes=tree$edge[idx,2]
  out=data.frame(node=nodes,start=nh[idx,1],end=nh[idx,2])
  if(length(nodes)>0){
    out$N=sapply(nodes,function(node) length(get_samples_in_clade(node,tree)))
  }else{
    out$N=integer()
  }
  out %>% filter(N>=min.count)
}
extract_timings=function(PD,threshold=50,min.count=2,tree.model="nb_tree",b.use.nodes=FALSE){
  mtree=PD$pdx$tree_ml
  treefitres=PD$fit[[tree.model]]$null
  if(!b.use.nodes){
    lineages=get_lineages(mtree,threshold=threshold,min.count=min.count )
  }else{
    lineages=get_lineages2(mtree,PD$nodes %>% filter(child_count>1) %>% (function(x) x$node))
  }
  if(dim(lineages)[1]<1){
    return(NULL)
  }
  lineages=lineages %>% left_join(PD$nodes[,c("node","driver3","driver")])
  do_tree_timings_raw(mtree,treefitres,lineages)
}

do_tree_timings_raw=function(mtree,treefitres,lineages){
  #lineages=get_lineages(mtree,threshold=threshold,min.count=min.count )
  tree=get_treefit_ci(treefitres)
  nh=nodeHeights(tree)
  idx=match(lineages$node,tree$edge[,2])
  ##Wild type = non-expanded 
  mt.clones=do.call("c",
                    lapply(lineages$node,
                           function(node) get_samples_in_clade(node,tree
                           )
                    )
  )
  wt.clones=setdiff(tree$tip.label,c("zeros",mt.clones))
  nh[,2]=round(nh[,2],6)
  ts=unique(nh[match(1:length(tree$tip.label),tree$edge[,2]),2],digits = 6)
  ts=ts[which(ts>1e-3)]
  cat("Found",length(ts),"timepoints")
  tmp=data.frame(samples=wt.clones,age=nh[match(match(wt.clones,tree$tip.label),tree$edge[,2]),2])
  nwtcolony=sapply(ts,function(t) length(which(tmp$age==t)))
  
  lineages$lower=tree$lower_lb95[idx]
  lineages$upper=tree$upper_ub95[idx]
  for(field in c("lower_lb95","lower_ub95","lower_median","upper_lb95","upper_median","upper_ub95")){
    lineages[[sprintf("t_%s",field)]]=tree[[field]][idx]
  }
  list(lineages=lineages)
}


get_treefit_ci=function(treefit){
  tree=treefit$ultratree
  tmp = tree
  post.dist=rstan::extract(treefit$fullres)
  nodetdist = matrix(apply(post.dist$ta, 1, function(x) {
    tmp$edge.length = x
    nodeHeights(tmp)[,1]
  }), nrow = length(tmp$edge.length))
  parent_bounds = apply(nodetdist, 1, function(x) quantile(x, prob = c(0.025, 
                                                              0.5, 0.975), na.rm = TRUE))
  nodetdist = matrix(apply(post.dist$ta, 1, function(x) {
    tmp$edge.length = x
    nodeHeights(tmp)[,2]
  }), nrow = length(tmp$edge.length))
  child_bounds = apply(nodetdist, 1, function(x) quantile(x, prob = c(0.025, 
                                                                       0.5, 0.975), na.rm = TRUE))
  
  tree$lower_lb95=parent_bounds[1,]
  tree$lower_median=parent_bounds[2,]
  tree$lower_ub95=parent_bounds[3,]
  tree$upper_lb95=child_bounds[1,]
  tree$upper_median=child_bounds[2,]
  tree$upper_ub95=child_bounds[3,]
  tree
}



collate_results=function(zz){
  cdf=getcolor_df(cohort %>% (function(x){x$patient=x$id;x}))
  allres=do.call("rbind",
                 lapply(names(zz),function(patient) {
                   out=zz[[patient]]$lineages
                   out$patient=patient
                   out}
                   )
                 )
  allres %>% rename(N="nchild") %>% left_join(cohort %>% dplyr::select(-patient),by=c("patient"="internal_id")) %>% left_join(cdf)
  
}

##First n terms of exponential as expressed as a proportion of the full exponential
expn=function(n,x){exp(-x)+sum(sapply(1:n,function(i) exp(-x+i*log(x)-lfactorial(i))))}

#p(t|m)~p(m|t)p(t) if we assume flat interval prior 0..To 
# then the cumulative density function is:
CDF=function(x,m=5,lambda=20,To=50) (sapply(x,function(x) (1-expn(m,lambda*x))/(1-expn(m,lambda*To))))
# So that we can straightforwardly include the development burst we can instead do the analysis in LAMBDA space
# and convert back to time at the end  (our flat prior is therefore now on 0-lambda*To)
CDF2=function(x,m=5,LAMBDA=50*20) (sapply(x,function(x) (1-expn(m,x))/(1-expn(m,LAMBDA))))
##In our standard development model we have excess mutations governed by the followin
logisticMean=function(L,k,midpoint,a,b){((L/k)*(log(1+exp(k*(b-midpoint)))-log(1+exp(k*(a-midpoint))))) / (b-a);}
totalLambda=function(lambda,t,t0=0){ (t-t0)*lambda+(t-t0)*logisticMean(149.2,-50.0,0.22,t0,t)}
##


simfromcdf=function(cdf,start,end,size){
  ##Simulates down to a granularity of 1/1e4 of the range.
  x=seq(start,end,length.out = 10000)
  cval=cdf(x)
  ##Find first value where cval>1e-
  idx=which(cval>1e-8)
  if(idx[1]>1){
    cval=cval[-(1:(idx[1]-1))]
    x=x[-(1:(idx[1]-1))]
  }
  rval=runif(size)
  idx=findInterval(rval,cval,rightmost.closed = TRUE,all.inside = TRUE)
  0.5*(x[idx]+x[idx+1])
}

sample_acq_time=function(m,lambda,To,N,t0=0){
  m=round(m)
  clambda=simfromcdf(function(x) CDF2(x,m=m,LAMBDA = lambda*(To-t0)),start = 0,end=lambda*(To-t0),size = N)
  sapply(clambda,function(clambda) { 
    #cat(clambda)
    uniroot(function(x) totalLambda(lambda,x,t0 = t0)-clambda,lower = max(1e-6,t0+1e-6),upper=To)$root}
  )
}

get_elapsed_time_tree_with_dev=function (tree, backgroundrate = NULL, 
                                         odf = 1) 
{
  N = length(tree$tip.label) + 1
  L = length(tree$edge.length)
  TT = max(tree$timestamp)
  idx.child = match(tree$edge[, 2], tree$edge[, 1])
  duration = ifelse(is.na(idx.child), TT - tree$tBirth, tree$tBirth[idx.child] - 
                      tree$tBirth)
  duration[which(tree$state == 0)] = 0
  t0=tree$tBirth/365
  t1=t0+duration/365+1e-6
  #browser()
  if (odf > 1) {
    tree$edge.length = rsimpop:::get_nb(n = L, meanmuts = totalLambda(backgroundrate,t1,t0), od = odf)
  }else {
    tree$edge.length = rpois(n=L, totalLambda(backgroundrate,t1,t0))
  }
  tree
}


## Mutation count 100 trees...
plot_all_trees=function(PDD,selres,mutcount){
  patients=unique(selres$patient)
  par(mfrow=c(length(patients),1),oma=c(5,5,5,5))
  
  for(patient in patients){
    PD=PDD[[patient]]
    tree=PD$pdx$tree_ml
    tmp=selres[which(selres$patient==patient),]
    tree$color=rep("grey",length(tree$edge.length))
    #if(dim(tmp)[1]>0){
      tree$color[match(tmp$node,tree$edge[,2])]=tmp$col
    #}else{
    #  next
    #}
    #tree$color=ifelse(tree$edge[,2] %in% selres$node[which(selres$patient==patient)],"black","grey")
    ultra=PD$fit$nb_tree$nullmodel$ultratree
    ultra$color=tree$color
    label=sprintf("%s:%s",PD$INTERNAL_ID,PD$OTHER_ID)
    plot_tree(tree,cex.label = 0,lwd=2);title(sprintf("%s:Mutation Burden Tree",label))
    #plot_tree(ultra,cex.label = 0,lwd=2);title(sprintf("%s:Time-Based Tree",label))
  }
  if(mutcount>0){
  mtext(sprintf("Expansions at %d mutations",mutcount),side = 3,cex=1.5,line = 2,outer = TRUE)
  }
}
## Based on ordering in meta-analysis - manual for now..
ogene=c("chr15_cnv","chr7_cnv","EIF6","RPL22","RPL5","TP53","PRPF8","GPR137B")


get_expansion_stats=function(tree){
  N=length(setdiff(tree$tip.label,"zeros"))
  sdi=sapply(1:200,function(i){
    lt=get_lineages(tree,i)
    n=sum(lt$N)
    p=c(rep(1/N,(N-n)),lt$N/N)
    p=p[p>0]
    sum(-p*log(p))
  })
  pct=sapply(1:200,function(i){
    lt=get_lineages(tree,i)
    n=sum(lt$N)
    n/N
  })
  list(N=N,sd_vs_count=sdi,pctclonal=100*pct)
}

report_lineages=function(PDD,threshold=75){
  cohort=data.frame(patient=PATIENTS,internal_id=sapply(PDD,function(PD) PD$INTERNAL_ID),supplier_id=sapply(PDD,function(PD) PD$OTHER_ID))
  cohort$NC=sapply(cohort$internal_id,function(x) length(PDD[[x]]$pdx$tree_ml$tip.label)-1)
  cohort$age_max=sapply(cohort$internal_id,function(x) max(PDD[[x]]$pdx$agedf$age_at_sample_exact[PDD[[x]]$pdx$agedf$age_at_sample_pcy>0.1]))
  cohort$age_min=sapply(cohort$internal_id,function(x) min(PDD[[x]]$pdx$agedf$age_at_sample_exact[PDD[[x]]$pdx$agedf$age_at_sample_pcy>0.1]))
  cohort=cohort %>% mutate(age_mean=0.5*(age_max+age_min)) ## unweighted mean
  lineages=do.call("rbind",lapply(PDD,function(PD){
    lng=get_lineages(PD$pdx$tree_ml,threshold = threshold)
    if(dim(lng)[1]==0){
      NULL
    }else{
      cbind(lng,internal_id=PD$INTERNAL_ID)
    }
  }))
  lineages %>% left_join(cohort)
}

get_driver_lineages=function(PD,threshold=75,mincount=2){
  l1=get_all_tree_drivers(PD$pdx,genes=GENES,cv=CV) %>% dplyr::select(-profile,-label2)
  l1=l1 %>% group_by(node) %>% summarise(label=paste(label,collapse=":"),cv=paste(cv,collapse=":"))
  if(dim(l1)[1]==0){
    l2=get_lineages(PD$pdx$tree_ml,threshold=threshold,min.count = mincount)
    if(dim(l2)[1]==0){
      return(list(driverlessclades=l2))
    }
    l2$overlap=0
    l2$internal_id=PD$INTERNAL_ID
    return(list(driverlessclades=l2))
  }else{
    tree=PD$pdx$tree_ml
    nh=nodeHeights(tree)
    idx=match(l1$node,tree$edge[,2])
    l1$start=nh[idx,1]
    l1$end=nh[idx,2]
    l1=cbind(l1,internal_id=PD$INTERNAL_ID,nchild=sapply(l1$node,function(x) length(get_samples_in_clade(x,PD$pdx$tree_ml))))
  }
  bycolony=do.call("rbind",lapply(l1$node,function(x) data.frame(node=x,colony=get_samples_in_clade(x,PD$pdx$tree_ml))))
  bycolony$internal_id=PD$INTERNAL_ID
  drivercolonies=do.call("c",lapply(l1$node,function(x) get_samples_in_clade(x,PD$pdx$tree_ml)))
  l2=get_lineages(PD$pdx$tree_ml,threshold=threshold,min.count = mincount)
  if(dim(l2)[1]>0){
    l2$overlap=sapply(l2$node,function(x) length(intersect(drivercolonies,get_samples_in_clade(x,PD$pdx$tree_ml))))
    l2$internal_id=PD$INTERNAL_ID
    l2=l2 %>% dplyr::rename("nchild"="N")
  }
  
  dat=list(driverclades=l1,driverlessclades=l2)#,bycolony=bycolony)
#  dat$driverlessclades2=dat$driverlessclades %>% 
#    left_join(dat$driverclades %>% dplyr::select(node,internal_id,label,nchild) %>% left_join(dat$bycolony,by=c("node","internal_id")),by=c("node","internal_id")) %>%
#    group_by(node,start,end,N,overlap,internal_id,label,nchild) %>% summarise(no=n())
  dat
}

combine_all_lineages=function(PDD,threshold=75,mincount=2){
  driverless=do.call("rbind",lapply(PDD,function(PD){
    dat=get_driver_lineages(PD,threshold = threshold,mincount = mincount)
    if(is.null(dat$driverlessclades)){
      return(NULL)
    }
    if(dim(dat$driverlessclades)[1]==0){
      return(NULL)
    }
    dat$driverlessclades
  }))
  #driverless=driverless %>% filter(overlap==0)
  drivers=do.call("rbind",lapply(PDD,function(PD){
    dat=get_driver_lineages(PD,threshold = threshold,mincount = mincount)
    if(is.null(dat$driverclades)){
      return(NULL)
    }
    if(dim(dat$driverclades)[1]==0){
      return(NULL)
    }
    dat$driverclades
  }))
  cohort=data.frame(patient=PATIENTS,internal_id=sapply(PDD,function(PD) PD$INTERNAL_ID))
  cohort$N=sapply(cohort$internal_id,function(x) length(PDD[[x]]$pdx$tree_ml$tip.label)-1)
  cohort$age_min=sapply(cohort$internal_id,function(x) min(PDD[[x]]$pdx$agedf$age_at_sample_exact[PDD[[x]]$pdx$agedf$age_at_sample_pcy>0.1]))
  cohort$age_max=sapply(cohort$internal_id,function(x) max(PDD[[x]]$pdx$agedf$age_at_sample_exact[PDD[[x]]$pdx$agedf$age_at_sample_pcy>0.1]))
  cohort$age_mean=0.5*(cohort$age_min+cohort$age_max)
  cohort=cohort %>% dplyr::select(-age_min,-age_max,-patient)
  rownames(drivers)=NULL
  rownames(driverless)=NULL
  drivers$status=ifelse(drivers$nchild>1,"driver_clonal_expansion","driver_singleton")
  driverless$status="nodriver_clonal_expansion"
  driverless$label="<NA>"
  driverless$cv="<NA>"
  drivers$overlap=NA
  ##browser()
  fields=c("internal_id","node","status","start","end","nchild","N","label","cv","age_mean","overlap")
  combo=rbind((drivers %>% left_join(cohort,by="internal_id"))[,fields],
              (driverless %>% left_join(cohort,by="internal_id"))[,fields])
  #list(drivers=drivers %>% left_join(cohort,by="internal_id"),driverless=driverless %>% left_join(cohort,by="internal_id"),c)
  combo
}

summarise_lineage_data=function(PDD,outstub="../export/lineage_summary"){
  s1=combine_all_lineages(PDD,threshold = 75) %>% dplyr::rename("patient"="internal_id")
  nodriveronly=s1 %>% filter(status=="nodriver_clonal_expansion")
  mixed=s1 %>% filter(is.na(overlap) | overlap==0)
  ## Driver based clonal expansions
  out=list(
    mixed=mixed,
  mixedsummary=mixed %>% group_by(patient,status,N,age_mean) %>% summarise(ncolony=sum(nchild),nevents=n()) %>% mutate(pct=ncolony/N),
  mixedsummarybypatient_all=mixed %>% group_by(patient,N,age_mean) %>% summarise(ncolony=sum(nchild),nevents=n()) %>% mutate(pct=ncolony/N),
  mixedsummarybypatient_ex_singleton=mixed %>% filter(status!="driver_singleton") %>% group_by(patient,N,age_mean) %>% summarise(ncolony=sum(nchild),nevents=n()) %>% mutate(pct=ncolony/N),
  nodriversummary=nodriveronly %>% group_by(patient,status,N,age_mean) %>% summarise(ncolony=sum(nchild),nevents=n()) %>% mutate(pct=ncolony/N)
  )
  out=lapply(out,function(x) x[order(x$age_mean),])
  labs=names(out)
  for(lab in labs){
    write.table(out[[lab]],sprintf("%s_%s.txt",outstub,lab),sep="\t",row.names = FALSE,quote=FALSE)
  }
  out
}

