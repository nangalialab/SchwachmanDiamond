
export_vcf=function(PD){
  format="REF_COUNT:ALT_COUNT:GENOTYPE"
  header=readLines("../data/VCFHEADER.txt")
  ## Add in samples
  samples=PD$pdx$meta$clones
  clones=PD$pdx$meta$clones_short
  sample_template="##SAMPLE=<ID=%s,Description=.,Accession=.,Platform=hiseq,Protocol=WGS,SampleName=%s,Source=.>"
  sampleLines=sprintf(sample_template,samples,samples)
  mtr=as.matrix(PD$pdx$dat$mtr[,clones])
  wtr=as.matrix(PD$pdx$dat$wtr[,clones])
  geno=as.matrix(PD$pdx$dat$geno_colony[,clones])
  geno=ifelse(is.na(geno),"?",as.character(geno))
  idx=1:length(mtr)
  entries=data.frame(matrix(sprintf("%d:%d:%s",wtr[idx],mtr[idx],geno),ncol=length(clones)),stringsAsFactors = FALSE)
  colnames(entries)=samples
  #info SOMATIC;MUT_TYPE=;GENE = ;PROTEIN_CHANGE=
  N=dim(PD$pdx$dat$details)[1]
  fix=function(x,default="-"){ifelse(is.na(x),"-",x)}
  INFO=with(PD$pdx$dat$details,sprintf("SOMATIC;MUT_TYPE=%s;GENE=%s;PROTEIN_CHANGE=%s;VC=%s;IS_IN_EXCLUDED_REGION=%d",TYPE,fix(GENE),fix(HGVS_PROTEIN),fix(VC),is_localx_excluded))
  info=with(PD$pdx$dat$details,data.frame(CHROM=Chrom,POS=Pos,ID=sprintf("%s_%d",PD$patient,1:N),REF=Ref,ALT=Alt,QUAL="99",FILTER="PASS",INFO=INFO,FORMAT=format,stringsAsFactors = FALSE))
  chrom=c(1:22,"X","Y")
  idx=order(match(info$CHROM,chrom),info$POS)
  out=cbind(info,entries)[idx,]
  rownames(out)=NULL
  headerLine=paste0("#",paste(colnames(out),collapse="\t"))
  c(header,sampleLines,headerLine,apply(out,1,paste,collapse="\t"))
}




report_variants_per_colony=function(PD){
  colonies=PD$pdx$meta$clones_short
  colonies_long=PD$pdx$cfg$LABEL[match(colonies,PD$pdx$cfg$SHORT_LABEL)]
  tree=PD$pdx$tree_ml
  details=PD$pdx$dat$details
  nodes=match(colonies,tree$tip.label)
  do.call("rbind",lapply(1:length(nodes),function(i){
    node=nodes[i]
    df=data.frame(node=get_parents(node,tree$edge))
    df$nshared=sapply(df$node,function(node) length(get_samples_in_clade(node,tree=tree)))
    vars=details %>% dplyr::select(Chrom,Pos,Ref,Alt,TYPE,GENE,CCDS,HGVS_PROTEIN,VC,is_localx_excluded,node) %>% inner_join(df,by="node")
    vars$colony=colonies_long[i]
    vars$internal_id=PD$INTERNAL_ID
    vars
  }))
}

