
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
do_tree_selection=function(PD,threshold=50,min.count=2,tree.model="nb_tree",b.use.nodes=FALSE){
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
  do_tree_selection_raw(mtree,treefitres,lineages)
}

do_tree_selection_raw=function(mtree,treefitres,lineages){
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
  sfit=lapply(1:length(lineages$node),function(i) {
    node=lineages$node[i]
    tmp=data.frame(samples=get_samples_in_clade(node,tree))
    tmp$age=nh[match(match(tmp$samples,tree$tip.label),tree$edge[,2]),2]
    nmutcolony=sapply(ts,function(t) length(which(tmp$age==t)))
    if(length(ts)>1){
      ##Fit using fit_Smt
      #browser()
      res=fit_Smt(nmutcolony = nmutcolony,nwtcolony = nwtcolony,ts=ts,
                  tmin=lineages$lower[i],tmax=lineages$upper[i],
                  LNmin = log(1e4),
                  LNmax=log(2e5),niter = 20000,
                  nchain = 3,
                  stan.control = list(adapt_delta=0.99))
    }else{
      res=fit_S(nmutcolony = length(tmp$samples),nwtcolony = length(wt.clones),T=ts,
                tmin=lineages$lower[i],tmax=lineages$upper[i],
                LNmin = log(1e4),
                LNmax=log(2e5),niter = 20000,
                nchain = 3,
                stan.control = list(adapt_delta=0.99))
      
    }
    ##Get S
    summary(res$res)$summary
  })
  list(lineages=cbind(lineages,t(sapply(sfit,function(x) x["S",c("2.5%","mean","97.5%","sd")])) %>% (function(x){colnames(x)=sprintf("S%s",colnames(x));x})),sfit=sfit)
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

# Estimating Selection

plot_selection_results=function(selres,mutcount,legpos="right",maxS=-1,b.add.id=FALSE){
 
  #tmp=tmp %>% left_join(cohort)
  selres[,grep("^S",colnames(selres))]=100*selres[,grep("^S",colnames(selres))]
  if(maxS<0){
  maxS=ceiling(max(selres$`S97.5%`))
  }
  par(mfcol=c(2,1))
  
  if(mutcount>0){
    extra=sprintf("from Expansions subsequent to %d Mutations",mutcount)
    selres=selres[order(selres$end),]
    tmp=unique(selres[,c("patient","col")])
    selres$pch=19
  }else{
    extra="for Identified Drivers"
    selres=selres[order(selres$driver3),]
    tmp=unique(selres[,c("driver3","col")])
    PCH=data.frame(patient=unique(selres$patient))
    PCH$pch=(0:25)[-4][1:length(PCH$patient)]
    selres=selres %>% left_join(PCH)
    #selres$pch=19
    b.add.id=TRUE
  }
  plot(NULL,xlim=c(-0.15*maxS,maxS*1.0),ylim=c(0,dim(selres)[1]+2),axes=F,xlab="Selection Coefficient (%)",ylab="",main=sprintf("Selection Coefficients %s",extra))
  #browser()
  segments(x0=selres$`S2.5%`,x1=selres$`S97.5%`,y0=1:dim(selres)[1],lwd=3,lend=2,col=selres$col)
  points(x=selres$Smean,y=1:dim(selres)[1],pch=selres$pch,col=selres$col,cex=1.5)
  if(!is.null(tmp$patient)){
    legend(legpos,tmp$patient,lwd=3,col=tmp$col,ncol=2,bty="n")
  }
  if(!is.null(tmp$driver3)){
    legend(legpos,tmp$driver3,lwd=3,col=tmp$col,ncol=2,bty="n")
  }
  text(x=-0.1*maxS,y=1:dim(selres)[1],labels=selres$N)
  text(x=-0.1*maxS,y=dim(selres)[1]+1.5,labels="Clade Size")
  if(b.add.id){
    #text(x=-0.2*maxS,y=1:dim(selres)[1],labels=selres$patient,pos = 4,offset = 0,cex=0.6)
    legend("topright",PCH$patient,pch=PCH$pch,ncol=length(PCH$pch),bty = "n")
  }
  abline(v=seq(0,maxS,50),lty="dotted")
  axis(side = 1,at=seq(0,maxS,50))
  plot(NULL,ylim=c(0,maxS),xlim=c(0,300),pch=19,ylab="Selection Coefficient (%)",xlab="Branch Timings (Mutation Count)")
  rect(xleft=0,xright = 50,ybottom = -1,ytop=maxS+1,col="pink",border=NA)
  points(selres$end,selres$Smean,pch=19)
  segments(x0=selres$start,x1=selres$end,y0=selres$Smean,lwd=0.5,col=selres$col)
  segments(x0=selres$end,y0=selres$`S2.5%`,y1=selres$`S97.5%`,lwd=2,col=selres$col)
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

plot_selection_vs_mutcount=function(selres,maxS=max(100*ceiling(selres$`S97.5%`))){
  plot(NULL,ylim=c(0,maxS),xlim=c(0,300),pch=19,ylab="Selection Coefficient (%)",xlab="Branch Timings (Mutation Count)")
  rect(xleft=0,xright = 50,ybottom = -1,ytop=maxS+1,col="pink",border=NA)
  points(selres$end,100*selres$Smean,pch=19,col=selres$colour)
  segments(x0=selres$start,x1=selres$end,y0=100*selres$Smean,lwd=0.5,col=selres$colour,lend=2)
  segments(x0=selres$end,y0=100*selres$`S2.5%`,y1=100*selres$`S97.5%`,lwd=3,col=selres$colour,lend=2)
  leg=unique(selres[order(selres$driver,decreasing = TRUE),c("driver","colour")])
  legend("right",leg$driver,col=leg$colour,lwd=3,cex=1)
}

plot_selection=function(selres,maxS=max(100*ceiling(selres$`S97.5%`))){
  tmp=selres #[order(selres$driver),]
  #maxS=max(tmp$`S97.5%`)
  #plot(NULL,xlim=c(0,100*ceiling(maxS)),ylim=c(0,dim(tmp)[1]+length(unique(tmp$driver))),yaxt="n",
  N=dim(tmp)[1]
  groups=unique(tmp$driver)
  plot(NULL,xlim=c(-maxS,maxS),ylim=c(0,N+length(groups)),yaxt="n",
       xlab="",
       ylab="",xaxt="n",cex.lab=1.5)
  axis(side = 1,at=seq(0,maxS,50))
  abline(v=seq(0,maxS,50),col="grey",lwd=0.5)
  #browser()
  mtext(text = "Selection Coefficient(%)",side=1,at = 0.5*maxS,cex=1,line=3)
  #abline(h=N+2*length(groups)+1,lty="dotted")
  k=N+length(groups)
  summarystat=data.frame(driver=groups,stringsAsFactors = FALSE)
  summarystat$mean=NA
  p=1
  for(group in groups){
    #abline(h=k+0.5,lty="dotted",col="grey")
    idx=which(tmp$driver==group)
    midy=mean(k-(1:length(idx))+1)
    for(i in idx){
      segments(x0=100*tmp$`S2.5%`[i],x1=100*tmp$`S97.5%`[i],y0=k,lwd=2,col=tmp$colour[i],lend=2)
      points(x=100*tmp$Smean[i],y=k,pch=18,cex=2,lwd=2)
      text(x=-1*maxS,y=k,sprintf("%-5s",tmp$internal_id[i]),cex=1,pos = 4)#tmp$driver[i]))
      text(x=-0.8*maxS,y=k,sprintf("%s",tmp$driver_orig[i]),cex=1,pos = 4)#tmp$driver[i]))
      
      k=k-1
    }
    #text(x=-0.95*maxS,y=midy,labels=group,col=tmp$colour[idx[1]],pos=4)
    
    if(length(idx)>1){
      fer=rma(yi = 100*tmp$Smean[idx],sei = 100*tmp$Ssd[idx])
    }else{
      row=tmp[idx,]
      fer=list(ci.lb=100*row$`S2.5%`,ci.ub=100*row$`S97.5%`,beta=100*row$Smean)
    }
    dwidth=0.2
    polygon(x=c(fer$ci.lb,fer$beta,fer$ci.ub,fer$beta),y=c(k,k+dwidth,k,k-dwidth),col=tmp$colour[idx[1]])
    text(x=-1*maxS,y=k,labels=group,col=tmp$colour[idx[1]],pos=4,font=2)
    abline(h=k+c(-0.5,0.5),col="black")#lty="dotted",col="grey")
    k=k-1
    #abline(h=k,lty="dotted",col="grey")
    #k=k-1
    summarystat$mean[p]=fer$beta
    p=p+1
  }
  #leg=unique(tmp[,c("driver","colour")])
  #legend("right",leg$driver,col=leg$colour,lwd=3)

  summarystat
  
}

mlog2=function(splus1){
  log(2)/log(splus1)
}
plot_selection_doublingtime=function(selres,maxS=ceiling(max(mlog2(1+selres$`S2.5%`)))){
  if(is.null(selres$flag)){
    selres$flag=FALSE
  }
  tmp=selres[order(selres$driver),]
  #maxS=max(tmp$`S97.5%`)
  #plot(NULL,xlim=c(0,100*ceiling(maxS)),ylim=c(0,dim(tmp)[1]+length(unique(tmp$driver))),yaxt="n",
  N=dim(tmp)[1]
  groups=unique(tmp$driver)
  plot(NULL,xlim=c(-maxS,maxS),ylim=c(0,N+2*length(groups)+1),yaxt="n",
       xlab="Doubling Time (Years)",
       ylab="",xaxt="n")
  axis(side = 1,at=seq(0,maxS,1))
  abline(v=seq(0,maxS,50),col="grey",lwd=0.5)
  #browser()
  
  #abline(h=N+2*length(groups)+1,lty="dotted")
  k=N+2*length(groups)
  summarystat=data.frame(driver=groups,stringsAsFactors = FALSE)
  for(group in groups){
    #abline(h=k+0.5,lty="dotted",col="grey")
    idx=which(tmp$driver==group)
    midy=mean(k-(1:length(idx))+1)
    for(i in idx){
      segments(x0=mlog2(1+tmp$`S2.5%`[i]),x1=mlog2(1+tmp$`S97.5%`[i]),y0=k,lwd=2,col=tmp$colour[i],lend=2)
      points(x=mlog2(1+tmp$Smean[i]),y=k,pch=18,cex=2,lwd=2)
      text(x=-0.95*maxS,y=k,sprintf("%-5s      %s",tmp$internal_id[i],tmp$driver_orig[i]),cex=1,pos = 4)#tmp$driver[i]))
      if(tmp$flag[i]){
        text(x=0.98*maxS,y=k,"*")
      }
      k=k-1
    }
    #text(x=-0.95*maxS,y=midy,labels=group,col=tmp$colour[idx[1]],pos=4)
    idx=which(tmp$driver==group & !tmp$flag)
    if(length(idx)>1){
      fer=rma(yi =tmp$Smean[idx],sei = tmp$Ssd[idx])
      fer=list(ci.lb=mlog2(1+fer$ci.ub),ci.ub=mlog2(1+fer$ci.lb),beta=mlog2(1+fer$beta))
    }else if(length(idx)==1){
      row=tmp[idx,]
      fer=list(ci.ub=mlog2(1+row$`S2.5%`),ci.lb=mlog2(1+row$`S97.5%`),beta=mlog2(1+row$Smean))
    }
    if(length(idx)>0){
    dwidth=0.2
    polygon(x=c(fer$ci.lb,fer$beta,fer$ci.ub,fer$beta),y=c(k,k+dwidth,k,k-dwidth),col=tmp$colour[idx[1]])
    text(x=-0.95*maxS,y=k,labels=group,col=tmp$colour[idx[1]],pos=4)
    abline(h=k+c(-0.5,0.5),col="black")#lty="dotted",col="grey")
    }else{
      abline(h=k+c(-0.5,0.5),col="black")#lty="dotted",col="grey")
    }
    k=k-1
    #abline(h=k,lty="dotted",col="grey")
    k=k-1
  }
  #leg=unique(tmp[,c("driver","colour")])
  #legend("right",leg$driver,col=leg$colour,lwd=3)
  
  
}


plot_selection_logscale=function(selres,maxS=max(ceiling(log(1+selres$`S97.5%`)))){
  mlog2=function(x){log(x)}
  if(is.null(selres$flag)){
    selres$flag=FALSE
  }
  tmp=selres#[order(selres$driver),]
  #maxS=max(tmp$`S97.5%`)
  #plot(NULL,xlim=c(0,100*ceiling(maxS)),ylim=c(0,dim(tmp)[1]+length(unique(tmp$driver))),yaxt="n",
  N=dim(tmp)[1]
  groups=unique(tmp$driver)
  plot(NULL,xlim=c(-maxS,maxS),ylim=c(0,N+2*length(groups)+1),yaxt="n",
       xlab="Selection Coefficient",
       ylab="",xaxt="n")
  lticks=1:floor(exp(maxS))
  ticks=log(lticks)
  #browser()
  #labs=ifelse(round((exp(ticks)-1)) %in% c(0,1,2,5,10),sprintf("%3.0f",100*(exp(ticks)-1)),"")
  #axis(side = 1,at=log(1:floor(exp(maxS))),labels = sprintf("%3.0f",100*((1:floor(exp(maxS)))-1)))
  axis(side = 1,at=ticks,labels = FALSE)
  rticks=c(0,1,5,10)
  mtext(side = 1,line = 1,at = log(rticks+1),text = sprintf("%3.0f",100*rticks))
  abline(v=ticks,col="grey",lwd=0.5)
  #browser()
  
  #abline(h=N+2*length(groups)+1,lty="dotted")
  k=N+2*length(groups)
  summarystat=data.frame(driver=groups,mean=NA,S_lb95=NA,S_ub95=NA,N=0,stringsAsFactors = FALSE)
  p=1

  for(group in groups){
    #abline(h=k+0.5,lty="dotted",col="grey")
    idx=which(tmp$driver==group)
    midy=mean(k-(1:length(idx))+1)
    maxb=-1
    for(i in idx){
      segments(x0=mlog2(1+tmp$`S2.5%`[i]),x1=mlog2(1+tmp$`S97.5%`[i]),y0=k,lwd=2,col=tmp$colour[i],lend=2)
      points(x=mlog2(1+tmp$Smean[i]),y=k,pch=18,cex=2,lwd=2)
      text(x=-0.95*maxS,y=k,sprintf("%-5s      %s",tmp$internal_id[i],tmp$driver_orig[i]),cex=1,pos = 4)#tmp$driver[i]))
      if(tmp$flag[i]){
        text(x=0.98*maxS,y=k,"*")
      }
      k=k-1
      if(mlog2(1+tmp$Smean[i])>maxb){
        maxb=mlog2(1+tmp$Smean[i])
      }
    }
    #text(x=-0.95*maxS,y=midy,labels=group,col=tmp$colour[idx[1]],pos=4)
    idx=which(tmp$driver==group & !tmp$flag)
    if(length(idx)>1){
      fer=rma(yi =tmp$Smean[idx],sei = tmp$Ssd[idx])
      fer=list(ci.lb=mlog2(1+fer$ci.lb),ci.ub=mlog2(1+fer$ci.ub),beta=mlog2(1+fer$beta))
    }else if(length(idx)==1){
      row=tmp[idx,]
      fer=list(ci.lb=mlog2(1+row$`S2.5%`),ci.ub=mlog2(1+row$`S97.5%`),beta=mlog2(1+row$Smean))
    }
    if(length(idx)>0){
      dwidth=0.2
      polygon(x=c(fer$ci.lb,fer$beta,fer$ci.ub,fer$beta),y=c(k,k+dwidth,k,k-dwidth),col=tmp$colour[idx[1]])
      text(x=-0.95*maxS,y=k,labels=group,col=tmp$colour[idx[1]],pos=4)
      abline(h=k+c(-0.5,0.5),col="black")#lty="dotted",col="grey")
    }else{
      fer=list(beta=maxb,ci.lb=NA,ci.ub=NA)
      abline(h=k+c(-0.5,0.5),col="black")#lty="dotted",col="grey")
    }
    k=k-1
    #abline(h=k,lty="dotted",col="grey")
    k=k-1
    summarystat$mean[p]=exp(fer$beta)-1
    summarystat$S_lb95[p]=exp(fer$ci.lb)-1
    summarystat$S_ub95[p]=exp(fer$ci.ub)-1
    summarystat$N[p]=length(idx)
    p=p+1
  }
  #leg=unique(tmp[,c("driver","colour")])
  #legend("right",leg$driver,col=leg$colour,lwd=3)
  summarystat
  
}

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

