library("colorspace")
do_tree_paper_plot=function(PD,b.add.lines=FALSE, cex.annot=0,b.extend.colours=TRUE,ff=0.15,prioritise.group=NULL,b.add.bs=FALSE,min.clone.count=-1,vspace.reserve=NA,force.count=NA,...){
  lm=ff/(1-ff)
  ds=get_driver_scheme()
  dff=get_all_tree_drivers(PD$pdx,genes=GENES,cv = CV)
  dff$group=gsub(":.+","",dff$label)
  dff$group=gsub("_[A-Z]$","",dff$group)
  dff=dff %>% left_join(ds,by=c("group"))
  dff=dff[order(dff$group,dff$cv),]
  print(dff)
  if(b.add.bs){
    PD=add_bs(PD,get_nexus(PD))
  }
  tree=PD$pdx$tree_ml
  ID=""#PD$INTERNAL_ID #sprintf("%s:%s",PD$INTERNAL_ID,PD$OTHER_ID)# PD$LONG_DESC
  if(dim(dff)[1]==0){
    if(is.na(vspace.reserve)){
      vspace.reserve=ifelse(dim(dff)[1]>=10,0.8,0.3)
    }
    if(is.na(force.count)){
      force.count=ifelse(dim(dff)[1]>=10,-1,8)
    }
    pdx=set_color_by_age(PD$pdx)
    tree=plot_basic_tree(pdx,ID,
                         style="vertical",
                         genes=GENES,
                         cex.label=0,
                         cex.annot.label = 1,b.show.just.group=TRUE,
                         cv = CV ,left.margin.prop = lm,cex.terminal.dots = 0.5,seed=123576,legpos = NULL,vspace.reserve=vspace.reserve,b.extend.colours = b.extend.colours,...
    );
    cord=tree$coords[match(1:length(tree$tip.label),tree$edge[,2]),]
    points(cord$b1,tree$top-cord$a1,col=tree$tip.color,cex=1,pch=19)
    legend("topleft",legend = sprintf("%3.1f",pdx$age.df$age),col=pdx$age.df$color,pch=19,title = "Age at Sample",bty = "n")
  }else{
    
    nodes=dff$node
    labels=dff$label
    hm=matrix("white",ncol=length(tree$tip.label),nrow=length(nodes))
    colnames(hm)=tree$tip.label
    rownames(hm)=labels;
    for(i in 1:length(nodes)){
      tips=get_samples_in_clade(nodes[i],tree);hm[i,]=adjustcolor(dff$colour[i],alpha.f = 0.15)
      hm[i,match(tips,tree$tip.label)]=dff$colour[i]
    }
    print(hm)
    #browser()
    if(min.clone.count>0){
      cols=c("black","black")
      lineages=get_lineages(tree,min.clone.count)
      hsc=rep("white",length(tree$tip.label))
      filler=hsc
      if(dim(lineages)[1]>0){
        inf=lapply(lineages$node,function(node) get_samples_in_clade(node,tree))
        k=1
        for(x in inf){
          hsc[match(x,tree$tip.label)]=cols[(k %% length(cols))+1]
          k=k+1
        }
      }
      hsc=rbind(hsc,filler)
      colnames(hsc)=colnames(hm)
      hm=rbind(hsc,hm)
      rownames(hm)[1:2]=c("Clonal Fraction","")
      txtcols=c("black","black",dff$colour)
    }else{
      txtcols=dff$colour
    } 
    pdx=set_color_by_age(PD$pdx)
    if(is.na(vspace.reserve)){
      vspace.reserve=ifelse(dim(dff)[1]>=10,0.8,0.3)
    }
    if(is.na(force.count)){
      force.count=ifelse(dim(dff)[1]>=10,-1,8)
    }
    tree=plot_basic_tree(pdx,ID,
                         style="vertical",
                         genes=GENES,
                         cex.label=0,
                         cex.annot.label = cex.annot,b.show.just.group=TRUE,
                         cv = CV ,left.margin.prop = lm,cex.terminal.dots = 1,seed=123576,legpos = NULL,
                         vspace.reserve=vspace.reserve,
                         b.extend.colours = b.extend.colours,
                         prioritise.group = prioritise.group,lwd=1.2,...
    );
    cord=tree$coords[match(1:length(tree$tip.label),tree$edge[,2]),]
    points(cord$b1,tree$top-cord$a1,col=tree$tip.color,cex=1,pch=19)
    legend("topleft",legend = sprintf("%3.1f",pdx$age.df$age),col=pdx$age.df$color,pch=19,title = "Age at Sample",bty = "n")
    tree=add_heatmap(tree,heatmap=hm,cex.label = 0.7,border=NA,txtcols=txtcols,force.count = force.count,b.add.lines = b.add.lines,font=1)
    #tree=add_heatmap(tree,heatmap=hm,cex.label = 0.7,border=NA,txtcols=txtcols,force.count = ifelse(dim(dff)[1]>10,-1,10),b.add.lines = b.add.lines,font=1)
    #lm=3*ff/(1-3*ff)
    if(b.add.bs){
      bs_labels2(tree,80)
    }
  }
  text(PD$INTERNAL_ID,x=length(tree$tip.label)/2,y=tree$top,pos = 3,cex=1.5,xpd=TRUE)
  tree
}



plot_all_trees_scaled=function(bw,mar=c(0,2,0,2),b.fixed.scale=TRUE,d=c(20,3,12,4),figure.path=sprintf("../export/figure_2%s.pdf",VERSION)){## box width fraction of total image height
  
  if(b.fixed.scale){
    m=c(600,300,650,950)
  }else{
    m=c(700,700,700,700)
  }
  ypm=(1-bw*sum(d))/sum(m)
  
  ## relative heights of rows
  h=m*ypm+bw*d
  h=h/sum(h)
  vsr=bw*d/(m*ypm)
  cat("vsr",vsr,"\n")
  pdf(figure.path,w=10,h=10,pointsize = 8)
  layout(matrix(c(1,1,1,1,1,1,1,1,
                  2,2,3,3,3,4,4,4,
                  5,5,6,6,6,7,7,7,
                  8,8,9,9,9,10,10,10
  ),ncol=8,byrow = TRUE),width=rep(0.125,8),heights = h)
  do_tree_paper_plot(PDD$SDS5,ff = 0.05,vspace.reserve = vsr[1],force.count=d[1],mar=mar)
  do_tree_paper_plot(PDD$SDS1,ff=0.18,vspace.reserve = vsr[2],force.count=3,ymax=300,mar=mar)
  do_tree_paper_plot(PDD$SDS2,ff=0.2,vspace.reserve = vsr[2],force.count=3,ymax=300,mar=mar)
  do_tree_paper_plot(PDD$SDS3,ff=0.2,vspace.reserve = vsr[2],force.count=3,ymax=300,mar=mar)
  do_tree_paper_plot(PDD$SDS4,ff=0.18,vspace.reserve = vsr[3],force.count=12,ymax=650,mar=mar)
  do_tree_paper_plot(PDD$SDS6,ff=0.2,vspace.reserve = vsr[3],force.count=12,ymax=650,mar=mar)
  do_tree_paper_plot(PDD$SDS7,ff=0.2,vspace.reserve = vsr[3],force.count=12,ymax=650,mar=mar)
  
  ## Change the TP53 colour scheme for biallelic vs mono-allelic.
  DS<<-get_driver_scheme()
  DS$colour[8]="#CAB2D6" ##lighten("#CAB2D6",0.2)
  get_driver_scheme_saved=get_driver_scheme
  get_driver_scheme<<-function() DS
  get_driver_scheme2<<-get_driver_scheme
  PD=PDD$SDS8
  do_tree_paper_plot(PD,ff=0.18,prioritise.group = "TP53",vspace.reserve = vsr[4],force.count=4,ymax=950,mar=mar)
  get_driver_scheme<<-get_driver_scheme_saved
  get_driver_scheme2<<-get_driver_scheme
  ##  Done fixing TP53 colour schem..
  do_tree_paper_plot(PDD$SDS9,ff=0.2,vspace.reserve = vsr[4],force.count=4,ymax=950,mar=mar)
  do_tree_paper_plot(PDD$SDS10,ff=0.2,vspace.reserve = vsr[4],force.count=4,ymax=950,mar=mar)
  dev.off()
}

plot_comp=function(b.abs.scale=FALSE){
  test=read.table("../export/lineage_summary_CH_mixed.txt",head=T,stringsAsFactors=FALSE)
  ## above doesn't include wild types so add them in..
  wt=test %>% group_by(patient,N) %>% summarise(ntot=sum(nchild)) %>% mutate(status="wild_type",nchild=N-ntot)
  ## Also some patients aren't include above because they are all wild-type add those back in
  patients=sprintf("SDS%d",1:10)
  missing.patients=patients[which(!patients %in% test$patient)]
  non.mut=data.frame(patient=missing.patients,status="wild_type")
  non.mut$nchild=sapply(PDD[non.mut$patient],function(x) length(x$pdx$meta$clones))
  non.mut$N=non.mut$nchild
  
  test2=rbind(test %>% dplyr::select(patient,nchild,status,N),wt %>% dplyr::select(patient,nchild,status,N),
              non.mut %>% dplyr::select(patient,nchild,status,N))
              ##data.frame(patient=c("SDS1","SDS3"),nchild=c(18,23),status="wild_type",N=c(18,23)))
  test2$patient=factor(test2$patient,levels=sprintf("SDS%d",1:10))
  
  pal=RColorBrewer::brewer.pal(9,"Paired")
  celltype=data.frame(status=c("driver_clonal_expansion","driver_singleton","nodriver_clonal_expansion","wild_type"),
                      longstatus=c("Clonal Expansion(Driver)","Singleton(Driver)","Clonal Expansion(No Driver)","Wild Type"),
                      col=c(pal[2],pal[1],"grey50","grey95")
  )
  test2=test2 %>% left_join(celltype)
  test2$longstatus=factor(test2$longstatus,rev(celltype$longstatus))
  sf=scale_fill_manual(breaks=celltype$longstatus,values=celltype$col)
  
  if(!b.abs.scale){
    p=ggplot(data=test2 %>% mutate(pct=nchild/N),aes(x=patient,y=100*pct,fill=longstatus)) +geom_bar(stat="identity",colour="black")+ylab("Cell Fraction(%)")+xlab("")+ guides(fill=guide_legend(title="")) + sf+ 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }else{
    p=ggplot(data=test2,aes(x=patient,y=nchild,fill=longstatus)) +geom_bar(stat="identity",colour="black")+ylab("Colony Count")+xlab("")+ guides(fill=guide_legend(title="")) + sf+ 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }
  p
}

do_acquisition_plot_SDS_merged=function(inf,
                                        extra="",xmax=50,b.add.text=FALSE,b.include.stem=FALSE,cex.text=0.5,b.add.bar=FALSE,i.order=1){
  inf2=inf 
  patients=unique(inf2$patient)
  ## Gestational Age at Birth
  ab=MEAN_AGE_AT_DELIVERY_DAYS/365.25
  ab2=(MEAN_AGE_AT_DELIVERY_DAYS+14)/365.25
  #manipulate ranges so that we can extend pre-birth period...
  prebirthprop=0.08
  inf=inf[order(i.order*inf$upper_median),]
  for(field in c("lower_median","lower_lb95","lower_ub95","upper_median","upper_lb95","upper_ub95","max_age_at_sample")){
    inf[[field]]=inf[[field]]-ab
    inf[[field]]=ifelse(inf[[field]]<0,xmax*prebirthprop*inf[[field]]/ab2,inf[[field]])
  }
  N=dim(inf)[1]
  par(mar=c(6,6,0,2)+0.1)
  plot(NULL,xlim=c(-xmax*prebirthprop,xmax*1.0),ylim=c(0,N+1),yaxt="n",xlab="",ylab="",bty="n",xaxt="n")
  axis(side = 1,at = seq(0,35,5),labels = seq(0,35,5))
  axis(side =1, at = -xmax*prebirthprop+xmax*prebirthprop*c(0,26/40,1),labels=c("0","26",""),cex.axis=0.6,las=1)
  axis(side =1, at = -xmax*prebirthprop+xmax*prebirthprop*c(13/40),labels=c("13"),cex.axis=0.6,las=1)
  mtext(side = 1,line = 4,text = "Age\n(Years)",at = xmax/2,cex=1)
  mtext(side = 1,line = 4,text = "Gestational Age\n(Weeks)",at = -prebirthprop*xmax,cex=1)
  
  w=0.1
  for(j in 1:N){
    if(b.add.bar){
      rect(xleft=0,xright=inf$max_age_at_sample[j],ybottom = j-0.5*0.9,ytop=j+0.5*0.9,col="gray95",border=NA)
    }else{
      shade_between(x = c(inf$upper_median[j],inf$max_age_at_sample[j]),
                    y1=c(j,j-0.5*inf$pct[j]),
                    y2=c(j,j+0.5*inf$pct[j]),
                    col=inf$col[j])
    }
  }
  
  rect(xleft=-prebirthprop*xmax,xright=0,ybottom=0,ytop=N+0.5-0.05,col="lightpink",border=NA,lwd=1)
  inf$pct=inf$child_count/inf$N
  if(b.include.stem){
    arrows(x0 =inf$lower_lb95,x1=inf$upper_ub95,y0=(1:N),col=inf$col,code = 3,angle=90,length = 0.05,lwd=2)
  }else{
    ww=0.05
    rect(xleft=inf$upper_lb95,xright=inf$upper_ub95,ybottom = (1:N)-ww,ytop=(1:N)+ww,col=inf$col,border=NA)
    segments(x0=inf$upper_lb95,y0=(1:N)-2*ww,y1=(1:N)+2*ww,col=inf$col,lwd=1)
    segments(x0=inf$upper_ub95,y0=(1:N)-2*ww,y1=(1:N)+2*ww,col=inf$col,lwd=1)
  }
  if(!b.add.text){
    tmp=inf %>% group_by(group,col) %>% summarise(N=n()) %>% (function(x) x[dim(x)[1]:1,])
    legend(x=27,y=N+1,gsub("_"," ",tmp$group),col=tmp$col,pch=15,cex=1.5,bty = "n",pt.cex=2)
    mtext(text = "Full Cohort",side = 2,line=1,las=0,cex=1)
  }
  return(inf2)
  
}

do_acquisition_plot_SDS_within_patient=function(inf,
                                                extra="",xmax=50,b.add.text=FALSE,b.include.stem=FALSE,cex.text=0.5,b.add.bar=FALSE,b.cluster.top=FALSE,b.add.axis=TRUE){
  
  inf=add_within_patient_layout(inf,b.cluster.top)
  inf2=inf 
  patients=unique(inf2$patient)
  N=length(patients)
  ## Gestational Age at Birth
  ab=MEAN_AGE_AT_DELIVERY_DAYS/365.25
  ab2=(MEAN_AGE_AT_DELIVERY_DAYS+14)/365.25
  #manipulate ranges so that we can extend pre-birth period...
  prebirthprop=0.08
  for(field in c("lower_median","lower_lb95","lower_ub95","upper_median","upper_lb95","upper_ub95","max_age_at_sample")){
    inf[[field]]=inf[[field]]-ab
    inf[[field]]=ifelse(inf[[field]]<0,xmax*prebirthprop*inf[[field]]/ab2,inf[[field]])
  }
  par(mar=c(3,6,0,2)+0.1)
  if(!b.add.axis){
    par(mar=c(0,6,0,2)+0.1)
  }
  plot(NULL,xlim=c(-xmax*prebirthprop,xmax*1.0),ylim=c(0,3*length(patients)),yaxt="n",xlab="",ylab="",bty="n",xaxt="n")
  if(b.add.axis){
    axis(side = 1,at = seq(0,35,5),labels = seq(0,35,5))
    axis(side =1, at = -xmax*prebirthprop+xmax*prebirthprop*c(0,26/40,1),labels=c("0","26",""),cex.axis=0.6,las=1)
    axis(side =1, at = -xmax*prebirthprop+xmax*prebirthprop*c(13/40),labels=c("13"),cex.axis=0.6,las=1)
  }
  k=1
  width=0.05#0.1
  xgap=0.3##0.5*(0.5-width)
  kk=0
  pch.sample=2
  pch.diagnosis=18
  yb=0
  for(patient in patients){
    idx=which(inf$patient==patient)
    if(length(idx)>1){
      gap=3/(length(idx)+1)
    }else{
      gap=1.5
    }
    yb=3*kk+gap
    rect(xleft=0,xright=inf$max_age_at_sample[idx[1]],ybottom=3*kk+xgap,ytop=3*(kk+1)-xgap,col="gray95",border=NA,lwd=1)
    rect(xleft=-prebirthprop*xmax,xright=0,ybottom=3*kk+xgap,ytop=3*(kk+1)-xgap,col="lightpink",border=NA,lwd=1)
    thistop=3*(kk+1)-xgap
    ww=3-2*xgap
    minb=3*kk+xgap
    if(b.add.bar){
      XR=38
      rect(xleft=35,xright=XR,ybottom=3*kk+xgap,ytop=3*(kk+1)-xgap,col="gray90",border=NA,lwd=1)
    }
    ww2=0.025
    for(j in idx){
      yb=3*kk+xgap+ww*inf$e1[j]
      if(b.include.stem){
        arrows(x0 =inf$lower_lb95[j],x1=inf$upper_ub95[j],y0=yb,col=inf$col[j],code = 3,angle=90,length = 0.05,lwd=2)
      }else{
        wi=0.05
        rect(xleft=inf$upper_lb95[j],xright=inf$upper_ub95[j],ybottom = yb-wi,ytop=yb+wi,col=inf$col[j],border=NA)
        segments(x0=inf$upper_lb95[j],y0=yb-2*wi,y1=yb+2*wi,col=inf$col[j],lwd=1)
        segments(x0=inf$upper_ub95[j],y0=yb-2*wi,y1=yb+2*wi,col=inf$col[j],lwd=1)
      }
      if(b.add.text){
        text(x=inf$max_age_at_sample[idx[1]]+5,y=yb,labels = inf$driver_full[j],pos = 4,cex=cex.text)
      }
      if(b.add.bar){
        rect(xleft=35,xright=XR,ybottom=3*kk+xgap+ww*inf$eb1[j],ytop=3*kk+xgap+ww*inf$et1[j],col=inf$col[j],border=NA,lwd=1)
      }
    }
    if(b.add.bar){
      segments(x0 = XR,x1=XR+0.1,y0=3*kk+xgap+ww*seq(0,1,0.1),lwd=0.5,col="black",lend=2)
      segments(x0 = XR,x1=XR+0.4,y0=3*kk+xgap+ww*seq(0,1,0.5),lwd=0.5,col="black",lend=2)
    }
    ymid=3*kk+1.5
    mtext(text = patient,side = 2,line=1,at=ymid,las=0,cex=1)
    kk=kk+1
  }
  
  if(!b.add.text){
    tmp=inf %>% group_by(group,col) %>% summarise(N=n())
    yloc=3*kk/4 +3
    xxl=0.9*xmax
    offset=0.08*xmax#5
    f=width*8
    for(i in 1:dim(tmp)[1]){
      rect(xleft=xxl,xright=xxl+offset,ybottom=yloc+f*i,ytop=yloc+f*i+2*width,col=tmp$col[i],border=NA)
      text(x=xxl+offset,y=yloc+f*i+width,labels = tmp$group[i],pos=4,offset=1)
    }
  }
  inf2
}