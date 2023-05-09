##This reflects the manual curation of driver nodes.
require("treemut")
add_driver_node_info=function(PD){
  #browser()
  loh=sapply(PD$pdx$meta$LOH,function(x) x$LABEL)
  idx=grep("delEIF6",loh)
  if(length(idx)>0){
    cat("editing LOH info for EIF6 SV\n")
    for(i in idx){
      PD$pdx$meta$LOH[[i]]$LABEL=PD$pdx$meta$LOH[[i]] %>% (function(x) sprintf("EIF6_del:%s:%3.1fMb-%3.1fMb",x$chr,x$start/1e6,x$end/1e6))
    }
  }
  if(FALSE){
    ##Not now adding PR137B:1-13transloc for annotation purposes..
  if(PD$pdx$meta$prefix=="PDHW1563"){
    ##Add a SV CNA.
    if(length(PD$pdx$meta$CNA)==0){
      PD$pdx$meta$CNA[[1]]=list(LABEL="GPR137B:1-13transloc",chr=1,start=0,end=1,ploidy=2
                              ,samples=c("bw","y","lo0017","lo0021","x","lo0013","be","lo0036","g","aw","lo0008","bh","cg","lo0052","s","af","cf","bc","lo0055","lo0031","w"))
      PD$inf$meta$CNA=PD$pdx$meta$CNA
    }else{
      stop("unexpected CNA for PDHW1563")
    }
  }
  }
  dff=get_all_tree_drivers(PD$pdx,genes=GENES,cv = CV)
  if(dim(dff)[1]>0){
    dff$status=-1
    nodes=do.call("rbind",lapply(unique(dff$node),function(node) {
      idx=which(dff$node==node)
      driver=paste(dff$label[idx],collapse=",")
      data.frame(node=node,driver=driver,status=1,stringsAsFactors = FALSE)
    }))
    nodes=nodes %>% filter(!is.na(node))
    nodes$child_count=sapply(nodes$node,function(node) length(get_samples_in_clade(node=node,tree=PD$pdx$tree_ml)))
    nodes=nodes[order(nodes$child_count,decreasing = TRUE),]
  }else{
    nodes=data.frame(node=integer(),driver=character(),status=integer(),child_count=integer(),stringsAsFactors = FALSE)
  }
 
  PD$inf$meta$CNA=PD$pdx$meta$CNA
  PD$inf$meta$LOH=PD$pdx$meta$LOH
  PD$nodes=nodes
  if(dim(PD$nodes)[1] !=dim(unique(PD$nodes))[1]){
    browser()
  }
  PD
}

