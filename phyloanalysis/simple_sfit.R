require("rstan")
require("rsimpop")
sfit.stan="
data{
real tmin;
real tmax;
real T;
real LNmin;
real LNmax;
int nmutcolony; //Could be number of mutant colonies or reads.
int nwtcolony; //Could be number of wild type colonies or reads.
}

parameters {
real<lower=0.01,upper=3> s; //instantaneous growth rate-> exp(s)-1 per Year.
real <lower=tmin,upper=tmax> t0; //acquisition time
real <lower=0> p;
real <lower=LNmin,upper=LNmax> LN; //ln(popsiz)
}



model {
real P;
s ~ uniform(0.01,3);
t0 ~ uniform(tmin,tmax); //
LN ~ uniform(LNmin,LNmax);
p ~ exponential(s/(1+s));//Asympototic form of the Birth-Death clone size (*exp(-s*(T-t0))) [See Tavere] 
P=p*exp(s*(T-t0)-LN);//Rescale to absolute mutant population size
//We assume here that the mutant population is additive to the usual stable HSC population.This is 
//also convenient for keeping the binomial probability <1.
nmutcolony ~ binomial(nmutcolony+nwtcolony,P/(1+P));
}

generated quantities {
real S;
S=exp(s)-1;
}
"
sfitmt.stan="
data{
real tmin;
real tmax;
real LNmin;
real LNmax;
int ntimepoint;
real multiplier[ntimepoint];  //Support a mixture of colony sampling and targeted follow up sampling..
int nmutcolony[ntimepoint]; //Could be number of mutant colonies or reads.
int nwtcolony[ntimepoint]; //Could be number of wild type colonies or reads.
real ts[ntimepoint];
}

parameters {
real<lower=0.01,upper=3> s; //instantaneous growth rate-> exp(s)-1 per Year.
real <lower=tmin,upper=tmax> t0; //acquisition time
real <lower=0> p;
real <lower=LNmin,upper=LNmax> LN; //ln(popsiz)
}



model {
real P;
s ~ uniform(0.01,3);
t0 ~ uniform(tmin,tmax); //
LN ~ uniform(LNmin,LNmax);
p ~ exponential(s/(1+s));//Asympototic form of the Birth-Death clone size (*exp(-s*(T-t0))) [See Tavere] 
for(i in 1:ntimepoint){
  P=p*exp(s*(ts[i]-t0)-LN);//Rescale to absolute mutant population size
  //We assume here that the mutant population is additive to the usual stable HSC population.This is 
  //also convenient for keeping the binomial probability <1.
  nmutcolony[i] ~ binomial(nmutcolony[i]+nwtcolony[i],multiplier[i]*P/(1+P));
}
}

generated quantities {
real S;
S=exp(s)-1;
}
"

SFITMT_STAN_MODEL=stan_model(model_code=sfitmt.stan,model_name = "sfitmt")
SFIT_STAN_MODEL=stan_model(model_code=sfit.stan,model_name = "sfit")


#' Fits a selection coefficient to summary expansion data.
#'
#' @param nmutcolony Number of sampled mutant colonies  (could also be aggregate mutant read count on shared clonal variants)
#' @param nwtcolony Number os sampled wild type colonies (could also be aggregate WT read count)
#' @param tmin  Minimum acquisition time
#' @param tmax  Maximum acquisition time
#' @param T Age at which colonies are sampled
#' @param LNmin Log scale minimum nHSC
#' @param LNmin Log scale maximum nHSC
#' @return list with STAN results (res) and Posterior (posterior)
fit_S=function(nmutcolony,nwtcolony,tmin,tmax,T,LNmin=log(1e4),LNmax=log(1e5),niter=1000,nchain=3,stan.control=list(adapt_delta=0.95)){
  dat=list()
  dat$nmutcolony=nmutcolony
  dat$nwtcolony=nwtcolony
  dat$LNmin=LNmin
  dat$LNmax=LNmax
  dat$T=T
  dat$tmin=tmin
  dat$tmax=tmax
  stanr=sampling(SFIT_STAN_MODEL,
                 data=dat,iter = niter,control=stan.control,chains = nchain,cores = nchain)
  list(posterior=rstan::extract(stanr),
       res=stanr,
       dat=dat
  )
}
#' Fits a selection coefficient to summary expansion data.
#'
#' @param nmutcolony Number of sampled mutant colonies  (could also be aggregate mutant read count on shared clonal variants)
#' @param nwtcolony Number os sampled wild type colonies (could also be aggregate WT read count)
#' @param tmin  Minimum acquisition time
#' @param tmax  Maximum acquisition time
#' @param T Age at which colonies are sampled
#' @param LNmin Log scale minimum nHSC
#' @param LNmin Log scale maximum nHSC
#' @return list with STAN results (res) and Posterior (posterior)
fit_Smt=function(nmutcolony,nwtcolony,ts,tmin,tmax,multiplier=rep(1.0,length(ts)),LNmin=log(1e4),LNmax=log(1e5),niter=1000,nchain=3,stan.control=list(adapt_delta=0.95)){
  dat=list()
  dat$nmutcolony=nmutcolony
  dat$nwtcolony=nwtcolony
  dat$ts=ts
  dat$multiplier=multiplier
  dat$ntimepoint=length(dat$ts)
  dat$LNmin=LNmin
  dat$LNmax=LNmax
  
  dat$tmin=tmin
  dat$tmax=tmax
  stanr=sampling(SFITMT_STAN_MODEL,
                 data=dat,iter = niter,control=stan.control,chains = nchain,cores = nchain)
  list(posterior=rstan::extract(stanr),
       res=stanr,
       dat=dat
  )
}

run_benchmark_sim=function(S=0.4,nyears=33){
  run_selection_sim(0.1,1/365,target_pop_size = 5e4,nyears_driver_acquisition = 10,nyears = nyears,fitness=log(1+S),minprop = 0.05)
}

#wrap rsimpop simulator inference
wrap_rsimpop_sfit=function(selsim,ncolony=50,LNmin=log(1e4),LNmax=log(1e5),b.use.perfect.priors=FALSE){
  st=get_subsampled_tree(selsim,ncolony)
  idx=which(st$events$driverid==1)
  if(length(idx)!=1){
    stop("No mutant clade")
  }
  node=st$events$node[idx]
  ts=st$events$ts[idx]
  idx=which(st$edge[,2]==node)
  tmin=st$tBirth[which(st$edge[,2]==node)]/365
  if(node>length(st$tip.label)+1){
    tmax=st$tBirth[which(st$edge[,1]==node)][1]/365
  }else{
    stop("Need at least 2 mutant clades to do inference..")
  }
  nmutcolony=length(get_samples_in_clade(node,st))
  nwtcolony=length(st$tip.label)-1-nmutcolony ##subtract 1 for the outgroup
  T=max(st$timestamp)/365
  Nhsc=selsim$cfg$compartment$popsize[2]
  if(b.use.perfect.priors){
    LNmin=log(0.99*Nhsc)
    LNmax=log(1.01*Nhsc)
    tmin=0.99*ts/365
    tmax=1.01*ts/365
  }
  zzz=fit_S(nmutcolony,nwtcolony,tmin,tmax,T,LNmin=LNmin,LNmax=LNmax,niter = 2000,nchain = 3,stan.control = list(adapt_delta=0.99))
  st=get_elapsed_time_tree(st)
  st$edge.length=st$edge.length/365
  list(S=quantile(zzz$posterior$S,prob=c(0.025,0.5,0.975)),sampledtree=st,dat=zzz$dat)
}

