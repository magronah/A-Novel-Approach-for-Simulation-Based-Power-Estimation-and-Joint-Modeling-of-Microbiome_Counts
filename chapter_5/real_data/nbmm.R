library(tidyr)
library(tidyverse)
library(tibble)
library(here)
library(DESeq2)
library(RhpcBLASctl)
library(glmmTMB)
library(NBZIMM)
library(AICcmodavg)
##############################################################
path   =   paste0(getwd(),"/real_data/")
source(paste0(path,"/fun.R"))
##############################################################
data        =   readRDS(paste0(path,"aut_data.rds"))
metadata    =   readRDS(paste0(path,"aut_metadata.rds"))
##############################################################
nam        =  "PRJNA644763"
count_dd   =   data[[nam]]
meta_dd    =   metadata[[nam]] %>% 
  setNames(c("subject", "group"))

dd_filt      =   filter_fun(count_dd,meta_dd,abund_thresh=10, 
                            sample_thresh=5)
###########################################################
pp   =  deseqfun(dd_filt,meta_dd,alpha_level=0.1,ref_name="ASD")
###########################################################
countdata  =  pp$data$countdata
meta_dd    =  pp$data$meta_data

meta_dd$dummy = factor(1)
meta_dd$normalizer =   sizeFactors(pp$object)
dd   =   data.frame(t(countdata))
###########################################################
nbmm_mod    =   mms(y = dd, fixed = ~group + offset(normalizer),
                   random =  ~ 1|dummy,
                   data = meta_dd, 
                   method = "nb",
                   niter = 100)
#saveRDS(mod_list, file=paste0(path,"autism/data/nbmm_runtime",nam,".rds"))
saveRDS(nbmm_mod, file=paste0(path,"autism/data/nbmm_mod",nam,".rds"))
###########################################################
par_nbmm = lapply(nbmm_mod$fit, function(x){
  list(beta = fixef(x), 
       sd   = log(as.numeric(VarCorr(x)["(Intercept)", "StdDev"])),
      theta = log(x$theta))
})
###########################################################
res   =  list()
num_params <- length(par_nbmm[[1]]) + 1 
n  =  names(nbmm_mod$fit)

for(i in 1:length(n)){
  y  =  dd[,n[i]]

  modA   =   glmmTMB(y ~ group + (1 | dummy) + offset(normalizer), 
                     data = meta_dd, 
                     family = nbinom2())
  
  
  modB1  =   glmmTMB(y ~ group + (1 | dummy) + offset(normalizer), 
                     data   = meta_dd, 
                     family = nbinom2(),
                     doFit  = FALSE)
  
  modB2 <- fitTMB(modB1, doOptim = FALSE)
  
  pars <- modB2$env$par
  f   =  par_nbmm[[i]]
  pars[names(pars) == "beta"] <- f$beta
  pars[names(pars) == "theta"] <- f$sd
  pars[names(pars) == "betadisp"] <- f$theta
  
  ## drop RE parameters
  pars <- pars[names(pars) != "b"]
  
  par2 <- with(modA$obj$env, last.par.best[-random])
  
  num_params =  attr(logLik(modA), "df")  
  res[[i]]  = list(nbmm_loglik    = modB2$fn(pars),
                   glmmTMB_loglik = modB2$fn(par2), 
                   params =  cbind(NBZIMM = pars, glmmTMB = par2),
                   AIC    = 2*(modB2$fn(pars)) +  2*(num_params),
                   mod    =   modA)
}
#####################################################################
diff = setdiff(colnames(dd), n)
res2 = list()
for(i in 1:length(diff)){
  y  =  dd[,diff[i]]
  
  mod   =   glmmTMB(y ~ group + (1 | dummy) + offset(normalizer), 
                     data = meta_dd, 
                     family = nbinom2())
  
  loglike = as.numeric(logLik(mod))
  
  res2[[i]] = list(mod = mod,
                   AIC    =  -2*(loglike) + 2*(num_params),
                   loglike = loglike)
}
####################################################################
nbmm_aic   = c(AIC  =   sum(unlist((lapply(res, function(x){x[["AIC"]]})))) +
                        sum(unlist((lapply(res2, function(x){x[["AIC"]]})))),
            loglik  =   -sum(unlist((lapply(res, function(x){x[["nbmm_loglik"]]})))) +
                        sum(unlist((lapply(res2, function(x){x[["loglike"]]}))))
)
saveRDS(nbmm_aic, file=paste0(path,"autism/data/nbmm_aic",nam,".rds"))
