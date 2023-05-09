args = commandArgs(trailingOnly=TRUE)#
if(length(args)<2){
	stop("Rscript plot_tree <patient prefix(e.g. PD6646d)> <pdf plot width> <do per colony|FALSE> <tree method(mpboot|scite_t|sifit> <remove.germline|TRUE>")
}
#source("trees.R")
source("load_and_annotate_tree.R")
source("trees_qc.R")
library("phytools")
#source("plot_tree_annots.R")
patient=args[1]
basewidth=as.integer(args[2])
ttype="mpboot"
b_per_colony=FALSE
b_remove_germline=TRUE
if(length(args)>2){
b_per_colony=as.logical(args[3])
}
if(lengths(args)>3){
  ttype=args[4]
  if( !(ttype %in% c("mpboot","scite_t","sifit"))){
    stop("invalid tree method")
  }
  if(length(args)>4){
  	b_remove_germline=as.logical(args[5])
  }
}
bs=system("ls ../results/*.splits.nex",intern = T)
bs=bs[grep(patient,bs)]
PD=get_phylo_dat(patient,ttype)
if(length(PD)>1){
	stop("Multiple runs of available")
}
outdir=sprintf("../tree_results/%s/%s",patient,ttype)
dir.create(outdir,recursive = TRUE)
bycolony=sprintf("%s/%s_by_colony.pdf",outdir,patient)
treeqc=sprintf("%s/%s_treeqc.pdf",outdir,patient)
rdat=sprintf("%s/%s_tree.Rdat",outdir,patient)
#badloci=sprintf("%s/%s_badloci.txt",outdir,patient)
#writeLines(PD$)

PD=get_phylo_dat(patient,ttype,PD)
PD$pdx$meta=PD$inf$meta

PD=add_bs(PD,bs)
PD$pdx$meta$b.remove.germline=b_remove_germline
pdf(treeqc,w=basewidth,h=12);
extra_qc_plot(PD);
plot_embryonic(PD$pdx,c(2,1))
par(mfcol=c(2,1))
plot_dodgy_variants(PD)
dev.off()
cat("Finished basic QC\n")
inf=PD$inf
pdx=PD$pdx
save(inf,pdx,file = rdat)
if(b_per_colony){
pdf(bycolony,w=basewidth,h=12);
  plot_all_vaf(PD$pdx,patient,min.mean.vaf=PD$pdx$meta$MINVAF);
dev.off()
}

cat("Finished..")

