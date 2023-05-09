#TELOMERE_FILE="../data/telomeres.2019_11_25.grouped_by_lcm_status.txt"
#BAITSET_PILEUP_SNV="../data/PD4781a_PD4781f_snp_vaf_vaf.tsv"
#BAITSET_PILEUP_INDEL="../data/PD4781a_PD4781f_indel_vaf_vaf.tsv"
#PATIENTS=c("PD42191","PD43293" )
PATIENTS=c("PD43294",  "PD43975",  "PD43295",  "PD43293",  "PDHW1563" ,"PD45887", "PD44706" , "PD42191"  ,"PD45889", "PD45886" )
##Reference colonies  : for cn analysis
#REFCOLONIES=list(PD7271="c012",PD5182="dl78",PD5163="d46",PD5847="142",PD5179="aw",PD6629="d16",
#                 PD5117="br",PD6646="70m",PD4781="d19",PD5147="c024",PD6634="lo0015",PD9478="cg")
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



add_agedf=function(PD){
  tree=PD$pdx$tree_ml
  sample=PD$pdx$cfg$LABEL[match(PD$pdx$tree_ml$tip.label,PD$pdx$cfg$SHORT_LABEL)]
  if(length(which(is.na(sample)==1))){
    sample[is.na(sample)]="zeros"
  }else{
    stop("unexpected lookup error!")
  }
  agedf=data.frame(tip.label=tree$tip.label,sample=sample,age_at_sample_exact=rep(1e-6,length(tree$tip.label)))
  agedf$age_at_sample_pcy=agedf$age_at_sample_exact+(MEAN_AGE_AT_DELIVERY_DAYS)/365.25
  PD$pdx$agedf=agedf
  PD
}
