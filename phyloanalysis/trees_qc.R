#GENES=readLines("GENES.txt")
#CV=readLines("CV.txt")

check_filters=function(dat,min.counts=0,max.count=1e6){
  res=list()
  filters=colnames(dat$inf$snv$details)[grep("^filter",colnames(dat$inf$snv$details))]
  clones=dat$inf$meta$clones_short
  for(filter in filters){
    idx=with(dat$inf$snv,which(details[[filter]]==1 & details$counts >= min.counts & details$counts <= max.count))
    if(length(idx)>10){
      info=assign_to_tree(dat$pdx$tree_ml,dat$inf$snv$mtr[idx,clones],dat$inf$snv$dep[idx,clones])#,error_rate=p.err,maxits=maxits)
      res[[filter]]=info$summary#assign_to_tree(dat$inf$snv$mtr[idx,clones],dat$inf$snv$dep[idx,clones],dat$pdx$df)
    }
  }
  res
}

get_overlap=function(inf,type){
  filters=colnames(inf[[type]]$details)[grep("^filter",colnames(inf[[type]]$details))]
  F=as.matrix(inf[[type]]$details[,filters])
  cnames=gsub("^filter_","",colnames(F))
  colnames(F)=cnames
  print(data.frame(filter=filters,NA_COUNT=colSums(is.na(F))))
  F[is.na(F)]=0
  count_matrix=t(F) %*% F
  total=length(which(rowSums(F)>0))
  prop_matrix=count_matrix/total
  list(count_matrix=count_matrix,prop_matrix=prop_matrix,total=total,pass=length(which(rowSums(F)==0)))
}

plot_table=function(tab,fmt="%s",title="Mutation Overlap",cex.headings=0.8){
  plot(NULL,xlim=c(0,1),ylim=c(0,1),axes = FALSE,xlab="",ylab="",main=title)
  R=dim(tab)[1]
  C=dim(tab)[2]
  dr=1/(R+2)
  dc=1/(C+2)
  sdr=sort(seq(dr,by = dr,length.out = R+1),decreasing = TRUE)
  sdc=seq(dc,by=dc,length.out = C+1)
  text(x=sdc[1],y=sdr[-1],labels=colnames(tab),cex=cex.headings)
  text(x=sdc[-1],y=sdr[1],labels=rownames(tab),cex=cex.headings)
  for(i in 1:R){
    font=rep(1,C)
    font[i]=2
    text(x=sdc[-1],y=sdr[i+1],labels=sprintf(fmt,tab[i,]),font=font)
  }
  arrows(x0=sdc+dc/2,y0=0,y1=1,lend=2,length=0)
  arrows(y0=sdr+dr/2,x0=0,x1=1,lend=2,length=0)
}



extra_qc_plot=function(dat){
  prefix=dat$inf$meta$prefix
  min.mean.vaf=dat$inf$meta$MINVAF
  par(mfrow=c(2,1),oma = c(0, 0, 5, 0))
  overlap.snv=get_overlap(dat$inf,"snv")
  overlap.indel=get_overlap(dat$inf,"indel")
  plot_table(100*overlap.snv$prop_matrix,fmt="%3.1f",title=sprintf("SNV:Pairwise Filter Overlap (%% Filtered Variants):Total Fail=%s/Pass=%s",overlap.snv$total,overlap.snv$pass),cex.headings = 0.8)
  plot_table(100*overlap.indel$prop_matrix,fmt="%3.1f",title=sprintf("Indels:Pairwise Filter Overlap (%% Filtered Variants):Total Fail=%s/Pass=%s",overlap.indel$total,overlap.indel$pass))
  mtext(prefix, outer = TRUE, cex = 1.5)
  par(mfcol=c(2,1))
  plot_vaf_by_branch(dat$pdx,prefix,b_pooled=TRUE)
  par(mfcol=c(1,1))
  tree=plot_basic_tree(dat$pdx,label = sprintf("%s: SNV+Indels",prefix))
  par(mfcol=c(2,1))
  dat=add_excluded_region_flag(dat)
  tmp=apply_adjustment_germline(dat)
  depth=colMeans(dat$pdx$dat$dep[,setdiff(dat$pdx$tree_ml$tip.label,"zeros")],na.rm=TRUE)
  tree=plot_basic_tree(tmp$pdx_snv,label = sprintf("%s: SNV With No Branch Length Adjustment",prefix),genes=GENES,cv=CV)
  tree=plot_basic_tree(tmp$pdx_snv_adj,label = sprintf("%s: SNV With Branch Length Adjustment",prefix),bars=depth,genes=GENES,cv=CV)
  ta1=tmp$pdx_snv_adj
  ta1$tree_ml$edge.length=ta1$tree_ml$el.snv.local.filtered
  tree=plot_basic_tree(ta1,label = sprintf("%s: SNV With No Branch Length Adjustment (local filters)",prefix),genes=GENES,cv=CV)
  ta1$tree_ml$edge.length=ta1$tree_ml$el.snv.local.filtered/ta1$tree_ml$sensitivity.snv.local.filtered
  tree=plot_basic_tree(ta1,label = sprintf("%s: SNV With No Branch Length Adjustment (local filters)",prefix),genes=GENES,cv=CV)

  ta1$tree_ml$edge.length=ta1$tree_ml$el.snv.global.filtered
  tree=plot_basic_tree(ta1,label = sprintf("%s: SNV With No Branch Length Adjustment (global filters)",prefix),genes=GENES,cv=CV)
  ta1$tree_ml$edge.length=ta1$tree_ml$el.snv.global.filtered/ta1$tree_ml$sensitivity.snv.global.filtered
  tree=plot_basic_tree(ta1,label = sprintf("%s: SNV With No Branch Length Adjustment (global filters)",prefix),genes=GENES,cv=CV)



  layout(matrix(1:2,ncol=1),heights=c(0.8,0.2))
  tmp=expand_short_branches(dat$pdx)
  tree=plot_basic_tree(tmp,label = sprintf("%s with node numbers (short branches extended to %s mutations)",prefix,tmp$minmuts),genes=GENES,cv=CV)
  node_labels(tree)
  par(mfcol=c(1,1))
  #plot_chromosome_annotated_tree(dat$pdx,sprintf("%s: Chromosome composition check",prefix))
  plot_chromosome_annotated_tree_pdx(dat$pdx,sprintf("%s: Chromosome composition check",prefix))
  plot_vaf_by_branch(dat$pdx,prefix,b_pooled = FALSE,b_horizontal = T,exclude.chr = "X")
  layout(matrix(1:2,ncol=1),heights=c(0.8,0.2))
  tree=plot_pooled_vaf_tree(dat$pdx,sprintf("%s [Pooled VAF]\nMin VAF=%3.2f",prefix,min.mean.vaf),min.mean.vaf = min.mean.vaf)
  #browser()
  if(!is.null(dat$pdx$tree_ml$BS)){
    minmuts=round(0.02*max(nodeHeights(dat$pdx$tree_ml)))
    layout(matrix(1:2,ncol=1),heights=c(0.8,0.2))
    tmp=dat$pdx
    tmp$tree_ml$color=ifelse(tmp$tree_ml$edge.length<minmuts,"grey","black")
    tmp$tree_ml$edge.length=ifelse(tmp$tree_ml$edge.length<minmuts,minmuts,tmp$tree_ml$edge.length)
    tree=plot_basic_tree(tmp,sprintf("%s: With Bootstrap Support (Grey branches extended to %s Mutations)",prefix,minmuts))
    bs_labels(tree)
  }

  chk=check_filters(dat,min.counts = 3)
  par(mfrow=c(3,3),oma = c(0, 0, 7, 0))
  ##browser()
  tmp=lapply(names(chk),function(x) hist(chk[[x]]$pval,breaks=seq(0,1,0.01),main=x,xlab="P-Value (Null hypothesis: Variant In Tree)") )
  mtext(sprintf("%s:Shared Variants (count>= 3)\n How well do discarded variants fit the tree?\nP-Values",prefix), outer = TRUE, cex = 1.5)
  tmp=lapply(names(chk),function(x) hist(chk[[x]]$p_else_where,breaks=seq(0,1,0.01),main=x,xlab="Probability Elsewhere in Tree") )
  mtext(sprintf("%s:Shared Variants (count>= 3)\n How well do discarded variants fit the tree?\nProbability Elsewhere in Tree",prefix), outer = TRUE, cex = 1.5)
  #browser()

  chk=check_filters(dat,min.counts = 1,max.count = 1)
  tmp=lapply(names(chk),function(x) hist(chk[[x]]$pval,breaks=seq(0,1,0.01),main=x,xlab="P-Value (Null hypothesis: Variant In Tree)") )
  mtext(sprintf("%s:Private Variants (count=1)\n How well do discarded variants fit the tree?\nP-Values",prefix), outer = TRUE, cex = 1.5)
  while(!par('page')) plot.new()
  tmp=lapply(names(chk),function(x) hist(chk[[x]]$p_else_where,breaks=seq(0,1,0.01),main=x,xlab="Probability Elsewhere in Tree") )
  mtext(sprintf("%s:Private Variants (count=1)\n How well do discarded variants fit the tree?\nProbability Elsewhere in Tree",prefix), outer = TRUE, cex = 1.5)
  while(!par('page')) plot.new()
  chk=check_filters(dat,min.counts = 2,max.count = 2)
  tmp=lapply(names(chk),function(x) hist(chk[[x]]$pval,breaks=seq(0,1,0.01),main=x,xlab="P-Value (Null hypothesis: Variant In Tree)") )
  while(!par('page')) plot.new()
  mtext(sprintf("%s:Weakly Shared Variants (count=2)\n How well do discarded variants fit the tree?\nP-Values",prefix), outer = TRUE, cex = 1.5)
  tmp=lapply(names(chk),function(x) hist(chk[[x]]$p_else_where,breaks=seq(0,1,0.01),main=x,xlab="Probability Elsewhere in Tree") )

  mtext(sprintf("%s:Weakly Shared (count=2)\n How well do discarded variants fit the tree?\nProbability Elsewhere in Tree",prefix), outer = TRUE, cex = 1.5)


}

check_haplotypes=function(dat){
  for(loh in c(dat$inf$meta$LOH,dat$inf$meta$CNA)){
    idx=with(dat,which(inf$snv$details$Chrom==loh$chr & inf$snv$details$Pos>loh$start & inf$snv$details$Pos<loh$end & inf$snv$details$G1000_AC>10 & inf$snv$details$BGLOD==1))
    cat("found",length(idx),"SNPS for ",loh$LABEL,".\n")
    if(length(idx)>5){
      samps=intersect(loh$samples,dat$pdx$meta$clones_short)
      if(length(samps)>10){
        pcor=cor(dat$inf$snv$vaf_d[idx,samps],use="pairwise.complete.obs")
        pheatmap(1-pcor,main=loh$LABEL)
      }else{
        if(length(samps)>1){
          vafs=lapply(samps,function(x) dat$inf$snv$vaf_d[idx,x])
          names(vafs)=samps
          pairs(vafs,cex=0.5,main=loh$LABEL)
        }
      }
    }
  }

}


plot_embryonic=function(PD,dims=c(3,3)){

  idx=which(PD$pdx$dat$details$BGLOD==1  )
  if(length(idx)==0){
    cat("No embryonic variants found...\n")
    return(NULL)
  }else{

  }
  labels=with(PD$pdx$dat$details[idx,],sprintf("Embryonic Variant?:%s:chr%s:%s",GENE,Chrom,Pos))
  if(length(idx)>9){
    cat("too many variants:",length(idx),"please provide plot dims for par(mfcol)\n")
  }
  par(mfcol=dims)
  PD$pdx=extend_branches(PD$pdx)
  tree=plot_tree(PD$pdx$tree_ml,b_do_not_plot=TRUE)
  lapply(1:length(idx),function(k){
    z=PD$inf$meta$clones_short[which(PD$pdx$dat$geno_colony[idx[k],]>0)];
    tree$tip.color=ifelse(tree$tip.label %in% z,"red","grey");
    plot_tree(tree,cex.label = 0.8);title(labels[k])})

}


plot_dodgy_variants=function(pdat,pvalue_threshold=1e-10){
  idx=which(pdat$pdx$summary$pval<1e-10)
  cat("Found",length(idx),"dodgy fits..\n")
  for(i in idx){tree=plot_tree(pdat$pdx$tree_ml,b_do_not_plot = TRUE)
  z=pdat$inf$meta$clones_short[which(pdat$pdx$dat$geno_colony[i,]>0)]
  tree$tip.color=ifelse(tree$tip.label %in% z,"red","grey")
  plot_tree(tree,cex.label = 0.8);title(sprintf("Bad fit p<1e-10:%s:%s",pdat$pdx$dat$details$Chrom[i],pdat$pdx$dat$details$Pos[i]))
  }
  idx
}

plot_vaf_by_branch=function(pdx,label,cex.label=1,b_pooled=T,vtype=c("SNV"),b_horizontal=FALSE,exclude.chr=c()){

  dd=get_vaf_info_by_branch(pdx,vtype,exclude.chr = exclude.chr)

  if(!b_pooled){
    idx.order=order(sapply(dd,function(x) x$ns))
    dd=dd[idx.order]
    ns=sapply(dd,function(x) x$ns)
    avafs=lapply(dd,function(x) x$vaf)
    if(b_horizontal){
      par(mar=c(5,8,4,2)+0.1)
      bnames=sprintf("%s:%s [%s]",ns,pdx$tree_ml$edge[idx.order,2],sapply(avafs,length))
      las=1
      xlab="VAF"
      ylab=""
    }else{
      bnames=sprintf("%s:%s\n[%s]",ns,pdx$tree_ml$edge[idx.order,2],sapply(avafs,length))
      las=1
      xlab="Edge"
      ylab="VAF"
    }

    idx.keep=which(ns>1)
    suppressWarnings(boxplot(avafs[idx.keep],names=bnames[idx.keep],notch=T,
                             main=sprintf("%s:VAF Distribution By Edge [Shared Branches]",label),
                             xlab=xlab,ylab=ylab,horizontal=b_horizontal,pars=list(cex.axis=cex.label),las=las))
    if(b_horizontal){
      abline(v=0.5)
    }else{
      abline(h=0.5)
    }
    return(NULL)
  }
  dcnt=unique(sapply(dd,function(x) x$ns))
  dcnt=sort(dcnt[dcnt>0])
  avafs=lapply(dcnt,function(d){idx=which(sapply(dd,function(x) x$ns)==d);do.call("c",lapply(dd[idx],function(x) x$vaf))})
  nv=sapply(avafs,length)

  bnames=sprintf("%s\n[%s]",dcnt,nv)
  suppressWarnings(boxplot(avafs,names=bnames,notch=T,
                           main=sprintf("%s:VAF Distribution By No. Colonies in Clade",label),
                           xlab="No of colonies in clade",ylab="VAF",pars=list(cex.axis=cex.label)))


  abline(h=0.5)

}

get_vaf_info_by_branch=function(pdx,vtype=c("SNV","INDEL"),exclude.chr=c()){
  lapply(pdx$tree_ml$edge[,2],function(i) {
    idx=get_idx_for_node(pdx$dat,node = i)
    idx=intersect(idx,which(pdx$dat$details$TYPE %in% vtype & pdx$dat$details$Chrom %in% setdiff(c(1:22,"X"),exclude.chr) ))
    samples=setdiff(get_samples_in_clade(node=i,pdx$tree_ml),"zeros")
    #browser()
    list(
      vaf=as.numeric(pdx$dat$mtr[idx,samples]/pdx$dat$dep[idx,samples]),
      ns=length(samples),
      samp=samples
    )

  }
  )

}

plot_embryonic=function(pdx,dims=c(3,3)){

  idx=which(pdx$dat$details$keep_embryonic==1  )
  if(length(idx)==0){
    cat("No embryonic variants found...\n")
    return(NULL)
  }else{

  }
  labels=with(pdx$dat$details[idx,],sprintf("Embryonic Variant?:%s:chr%s:%s",GENE,Chrom,Pos))
  if(length(idx)>9){
    cat("too many variants:",length(idx),"please provide plot dims for par(mfcol)\n")
  }
  par(mfcol=dims)
  pdx=expand_short_branches(pdx)
  tree=plot_tree(pdx$tree_ml,b_do_not_plot=TRUE)
  lapply(1:length(idx),function(k){
    z=pdx$meta$clones_short[which(pdx$dat$geno_colony[idx[k],]>0)];
    tree$tip.color=ifelse(tree$tip.label %in% z,"red","grey");
    plot_tree(tree,cex.label = 0.8);title(labels[k])})

}

plot_dodgy_variants=function(pdat,pvalue_threshold=1e-10){
  idx=which(pdat$pdx$summary$pval<1e-10)
  cat("Found",length(idx),"dodgy fits..\n")
  for(i in idx){tree=plot_tree(pdat$pdx$tree_ml,b_do_not_plot = TRUE)
  z=pdat$inf$meta$clones_short[which(pdat$pdx$dat$geno_colony[i,]>0)]
  tree$tip.color=ifelse(tree$tip.label %in% z,"red","grey")
  plot_tree(tree,cex.label = 0.8);title(sprintf("Bad fit p<1e-10:%s:%s",pdat$pdx$dat$details$Chrom[i],pdat$pdx$dat$details$Pos[i]))
  }
  idx
}

add_bs=function(dat,nexus.file){
  nexus=read_nexus(nexus.file)
  bs=data.frame(BS=sapply(nexus$splits,function(x) x$BS),
                SAMPLE=sapply(nexus$splits,function(x) paste(sort(x$taxa),collapse=":")))

  clades=get_clade_keys(dat$pdx)
  idx1=match(clades$c1,bs$SAMPLE)
  idx2=match(clades$c2,bs$SAMPLE)
  idx=ifelse(is.na(idx1),idx2,idx1)
  if(length(which(is.na(idx)))>0){
    cat("can't find all supporting bipartition setting BS to 0")
    # print(clades[idx,])
  }
  dat$pdx$tree_ml$BS=bs$BS[idx]
  dat$pdx$tree_ml$BS=ifelse(is.na(dat$pdx$tree_ml$BS),0,dat$pdx$tree_ml$BS)
  dat
}

bs_labels=function(tree){
  idx=which(tree$edge[,2]>length(tree$tip.label))
  text(tree$BS[idx],x=tree$coords$b1[idx],y=tree$top-tree$coords$a1[idx],col=ifelse(tree$BS[idx]<80,"red","blue"))
}

read_nexus=function(nexus.file){
  lines=readLines(nexus.file)
  idx1=which(lines=="TAXLABELS")+1
  idx2=which(lines=="END; [Taxa]")-2
  rng=idx1[1]:idx2[1]
  samples=gsub("'","",gsub("^\\[[0-9]+\\] ","",lines[rng]))
  idx1=which(lines=="MATRIX")+1
  idx2=which(lines=="END; [Splits]")-2
  rng=idx1[1]:idx2[1]
  splits=strsplit(lines[rng],"\t")
  #browser()
  splits=lapply(splits,function(x) {
    BS=as.numeric(x[2])
    taxa=unlist(strsplit(gsub(",$","",x[3])," "))
    taxa=taxa[which(nchar(taxa)>0)];list(BS=BS,taxa=samples[as.integer(taxa)])})
  list(taxa=samples,splits=splits)
}

get_clade_keys=function(pdx){
  tips=pdx$tree_ml$tip.label
  key=t(sapply(pdx$tree_ml$edge[,2],function(node){
    s=get_samples_in_clade(node = node,tree = pdx$tree_ml)
    c(c1=paste(sort(s),collapse=":"),c2=paste(sort(setdiff(tips,s)),collapse=":"))
  }))
  data.frame(key)
}

get_all_excluded=function(PDD=NULL,fname="loh_cna.exclude.txt",b_use_cache=TRUE){
  if(file.exists(fname) && b_use_cache){
    cat("Reading pre-generated exclude file",fname,"\n")
    return(read.table(fname,sep="\t",head=TRUE,colClasses = c("character","integer","integer"),stringsAsFactors = FALSE))
  }
  if(is.null(PDD)){
    stop(sprintf("need to specify PDD to generate %s",fname))
  }
  allf=do.call("rbind",lapply(PDD,function(x) get_cna_loh(x$inf)))
  #Exclude X & Y

  allf=rbind(allf,data.frame(chr=c("X","Y"),start=c(0,0),end=c(3e8,3e8)))
  grinds=GRanges(allf$chr, IRanges(start = allf$start, end=allf$end))
  chrinfo=get_chr_info()
  allgrinds=GRanges(chrinfo$chr, IRanges(start = chrinfo$start, end=chrinfo$end))
  df=as.data.frame(intersect(union(grinds,grinds),allgrinds))#Must be a more elegant way of merging intervals - but this works
  colnames(df)[1]="chr"
  df=df[,c("chr","start","end")]
  ##Now mark these
  if(b_use_cache){
    write.table(df,file=fname,sep="\t",quote=FALSE,row.names = FALSE)
  }
  df
}



