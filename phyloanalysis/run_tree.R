argsx = commandArgs(trailingOnly=TRUE)#
if(length(argsx)<4){
	stop("Rscript runTree.R <patient prefix(e.g. PD6646d)> <type(scite_t|sifit)> <label> <niter> <runs> <mem.gb> <queue> <input_dir>")
}

source("trees.R")
source("export.R")
args=argsx
runs=10
mem.gb=2
queue="normal"
qc.dir=NULL
patient=args[1]
stype=args[2]
label=args[3]
niter=as.numeric(args[4])
if(length(args)>4){
runs=as.integer(args[5])
}
if(length(args)>5){
mem.gb=as.integer(args[6])
}
if(length(args)>6){
queue=args[7]
}
if(length(args)>7){
qc.dir=args[8]
}
dir.create("../results",showWarnings=FALSE)
res=infer_tree(patient,type = stype,label=label,niter=niter,runs = runs,mem.gb = mem.gb,queue=queue,qc.dir=qc.dir)
save(res,file=sprintf("%s_%s_%s.Rdat",res$stub,niter,runs))
rdat=sprintf("%s_%s_%s.Rdat",res$stub,niter,runs)
dat=get_dat(patient,qc.dir=qc.dir)
dat=glue_infered_tree(rdat,pdat = dat)
save(dat,file=sprintf("%s_%s_%s_info.Rdat",res$stub,niter,runs))

