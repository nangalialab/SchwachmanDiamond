remap_details=function(details){
  dd=details
  colnames(dd)[match(c("Chrom","Pos","Ref","Alt","GENE","VC"),colnames(dd))]=c("chr","position","ref","mut","gene","vc")
  dd
}
##Add LOH
add_loh=function(inf,chr,start,end,samples){
  flag=inf$details$chr==chr & inf$details$position>=start & inf$details$position<=end
  ##Exclude potential germline variants that are now homozygous reference.
  idx.sample=match(samples,colnames(inf$gt))
  #browser()
  for(i in idx.sample){
    inf$gt[flag & inf$gt[,i]==0,i]=NA
    inf$gt[flag & inf$gt[,i]==2,i]=1
  }
  ##Keep 
  if(is.null(inf$cfg$LOH)){
    inf$cfg$LOH=list()
  }
  n=length(inf$cfg$LOH)
  inf$cfg$LOH[[n+1]]=list(chr=chr,start=start,end=end,samples=samples,label=sprintf("LOH_%d",n+1))
  inf
}

add_key=function(inf){
  if(!is.null(inf$details$key)){
    return(inf)
  }
  out=NULL;
  cnt=data.frame(gene=c("x"),count=c(0));
  #gene=gsub(pattern = "\\.",replacement = "xgg",inf$details$gene)
  ##Replace - with _
  ##Replace . with _
  ##Replace leading digit with _<digit?
  gene=gsub("([a-zA-Z0-9])-","\\1_",inf$details$gene)
  gene=gsub("\\.","_",gene)
  gene=gsub("^([0-9])","_\\1",gene)
  
  for(g in gene){
    idx=which(cnt$gene==g)
    if(length(idx)>0){
      out=c(out,cnt$count[idx]+1)
      cnt$count[idx]=cnt$count[idx]+1
    }else{
      out=c(out,1)
      cnt=rbind(cnt,data.frame(gene=g,count=c(1)))
    }
  }
  nr=nchar(inf$details$ref)
  nm=nchar(inf$details$mut)
  code=sprintf("%s_%s",inf$details$ref,inf$details$mut)
  for(i in which(nr>1 | nm>1)){
    if(nr[i]>1 && nm[i]>1){
      code[i]=sprintf("indel_d%s_i%s",nr[i]-1,nm[i]-1)
    }else{
      if(nr[i]>1){
        code[i]=sprintf("del%s",nr[i]-1)
      }else{
        code[i]=sprintf("ins%s",nm[i]-1)
      }
      
    }
    
  }
  
  inf$details$key=ifelse(inf$details$gene==".",sprintf("chr%s_%s_%s",inf$details$chr,inf$details$position,code),sprintf("%s_%s",gene,out))
  if(length(unique(inf$details$key))!=length(inf$details$key)){
    tmp=as.data.frame(table(inf$details$key),stringsAsFactors=FALSE)
    dups=tmp[which(tmp$Freq>1),1]
    for(dup in dups){
      idx=which(inf$details$key==dup)
      inf$details$key[idx]=sprintf("%s_u%s",dup,1:length(idx))
    }
    if(length(inf$details$key)!=length(inf$details$key)){
      stop("Duplicates in key...")
    }else{
      warning("Fixed duplicates")
    }
    
  }
  rownames(inf$gt)=inf$details$key
  inf
}

get_profile=function(gtm){
  sprintf("P%s",apply(gtm,1,function(x){paste(ifelse(is.na(x),".",as.character(x)),collapse="")}))
}
get_gpeq=function(gpe,quality.threshold=0.01){
  sprintf("Q%s",apply(gpe,1,function(x){paste(ifelse(x>quality.threshold,0,1),collapse="")}))
}
get_gpeq9=function(gpe){
  sprintf("Q%s",apply(gpe,1,function(x){qq=floor(-log10(x));paste(ifelse(qq>9,9,qq),collapse="")}))
}
##These functions use a legacy data structure.
export=function(inf,out.stub){
  if(is.null(inf$gt)){
    stop("Need to finalise genotypes before exporting")
  }
  inf=add_key(inf)
  #browser()
  gtm=inf$gt
  inf$details$profile=get_profile(gtm)
  gtm[which(is.na(gtm))]=3
  inf$details$gpq=get_gpeq9(inf$gpe)
  details=inf$details
  gpe=inf$gpe
  m=dim(gtm)[2]
  for(loh in inf$cfg$LOH){
    geno=rep(0,m)
    geno[match(loh$samples,colnames(gtm))]=1
    gtm=rbind(gtm,geno)
    rownames(gtm)[dim(gtm)[1]]=loh$label
    gpe=rbind(gpe,rep(1e-10,m))
    detl=tail(details,1)
    detl[1,]=NA
    detl[,"chr"]=loh$chr
    detl[,"position"]=loh$start
    detl[,"profile"]=sprintf("P%s",paste(geno,collapse=""))
    if(!is.na(match(loh$label,details$key))){
      stop("LOH:Attempt to add duplicate key")
    }
    detl[,"key"]=loh$label
    detl[,"vc"]="loh"
    details=rbind(details,detl)
  }
  
  #browser()
  write.table(gtm,file =sprintf("%s_geno.txt",out.stub) ,sep = "\t",row.names = FALSE,col.names = FALSE)
  write.table(floor(-10*log10(gpe)),file =sprintf("%s_gpe_phred.txt",out.stub) ,sep = "\t",row.names = FALSE,col.names = FALSE)
  write.table(details,file =sprintf("%s_details.txt",out.stub) ,sep = "\t",quote = FALSE,row.names = FALSE)
  writeLines(colnames(gtm),con =sprintf("%s_sample_names.txt",out.stub))
  writeLines(paste(colnames(gtm),collapse=" "),con=sprintf("%s_sample_names_sifit.txt",out.stub))
  writeLines(rownames(gtm),con =sprintf("%s_gene_names.txt",out.stub))
  ##The following seems to be required to get MPNboot to not crash (binary characters don't work)
  vals=apply(gtm,2,function(x) paste(ifelse(x==3,"?",ifelse(x==0,"T","A")),collapse=""))
  fasta=sapply(1:dim(gtm)[2],function(i) sprintf(">%s\n%s",colnames(gtm)[i],vals[i]))
  writeLines(fasta,con =sprintf("%s.fa",out.stub))
}



