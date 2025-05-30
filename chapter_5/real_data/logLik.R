library(MASS)
library(NBZIMM)
library(DESeq2)
library(tibble)
library(dplyr)
path        =   paste0(getwd(),"/reproducible/new_sim/real_data/")
##############################################################
source(paste0(path,"fun.R"))
source("load_glmmTMB.R")
library(glmmTMB)
source(paste0(path,"autism/load_data.R"))

rr_mod      =   readRDS(paste0(path,"autism/data/rr_mod.rds"))
us_mod      =   readRDS(paste0(path,"autism/data/us_mod.rds"))
nbmm_mod    =   readRDS(paste0(path,"autism/data/nbmm_mod.rds"))
deseq_mod    =   readRDS(paste0(path,"autism/data/deseq_mod.rds"))
zinbmm_mod  =   readRDS(paste0(path,"autism/data/zinbmm_mod.rds"))
##############################################
par_list = lapply(nbmm_mod$fit, function(x){
  pars  =   list(beta = fixef(x), 
                 sd   = log(as.numeric(VarCorr(x)["(Intercept)", "StdDev"])),
                 theta = log(x$theta))
})

dd = t(dd_filt)
met_dd$dummy  =  factor(1)
res =  list()

num_params <- length(par_list[[1]]) + 1  # +1 for the dispersion parameter


for(i in 1:length(nbmm_mod$fit)){
  
  n  =  names(nbmm_mod$fit)
  y  =  dd[,n[i]]
  
  ## fit model in glmmTMB
  modA   =   glmmTMB(y ~ group + (1 | dummy), 
                     data = met_dd, 
                     family = nbinom2())
  
  
  ## set-up model in glmmTMB
  modB1  =   glmmTMB(y ~ group + (1 | dummy), 
                     data = met_dd, 
                     family = nbinom2(),
                     doFit = FALSE)
  
  modB2 <- fitTMB(modB1, doOptim = FALSE)
  
  # names(modB2$env$par)
  pars <- modB2$env$par
  f   =  par_list[[i]]
  pars[names(pars) == "beta"] <- f$beta
  pars[names(pars) == "theta"] <- f$sd
  pars[names(pars) == "betadisp"] <- f$theta
  
  
  ## drop RE parameters
  pars <- pars[names(pars) != "b"]
  
  par2 <- with(modA$obj$env, last.par.best[-random])
  
  
  res[[i]]  = list(nbmm_loglik   = modB2$fn(pars),
       glmmTMB_loglik = modB2$fn(par2), 
       params =  cbind(NBZIMM = pars, glmmTMB = par2),
       AIC    =  2*(num_params - (-modB2$fn(pars)) ))
}

sum_a <- sum(sapply(res, function(x) x$AIC))
AIC(rr_mod)
logLik(rr_mod)
# form      =   count ~ 1 + us(1 + group|taxon) +  
#   rr(0 + taxon | subject,2) + 
#   offset(normalizer)


# rank   =   2
# ntaxa   =  max(dim(dd_filt))
# nsubj   =  min(dim(dd_filt))
# 
# num_params2 <- length(fixef(rr_mod)$cond)  +  
#                       3 * ntaxa            +  #us(1 + group | taxon): 3 * ntaxa
#                  ntaxa*rank  - choose(rank,2) + 
#                1
#               
# logLik(rr_mod)
# fixef(rr_mod) + 1
# AIC(rr_mod)


par_list = lapply(zinbmm_mod$fit, function(x){
  pars  =   list(beta = fixef(x), 
                 sd   = log(as.numeric(VarCorr(x)["(Intercept)", "StdDev"])),
                 theta = log(x$theta))
})


num_params <- length(par_list[[1]]) + 
              length(zinbmm_mod$fit[[1]]$zi.fit$coefficients) +
              1      # +1 for the dispersion parameter

for(i in 1:length(zinbmm_mod$fit)){
i=1  
  n  =  names(zinbmm_mod$fit)
  y  =  dd[,n[i]]
  
  ## fit model in glmmTMB
  modA   =   glmmTMB(y ~ group + (1 | dummy), 
                     ziformula = ~1,
                     data = met_dd, 
                     family = nbinom2())
  
  
  ## set-up model in glmmTMB
  modB1  =   glmmTMB(y ~ group + (1 | dummy), 
                     ziformula = ~1,
                     data = met_dd, 
                     family = nbinom2(),
                     doFit = FALSE)
  
  modB2 <- fitTMB(modB1, doOptim = FALSE)
  
  names(modB2$env$par)
  pars <- modB2$env$par
  f   =  par_list[[i]]
  pars[names(pars) == "beta"] <- f$beta
  pars[names(pars) == "theta"] <- f$sd
  pars[names(pars) == "betadisp"] <- f$theta
  pars[names(pars) == "betazi"] <- 
    as.numeric(zinbmm_mod$fit[[i]]$zi.fit$coefficients)
  
  
  ## drop RE parameters
  pars <- pars[names(pars) != "b"]
  
  par2 <- with(modA$obj$env, last.par.best[-random])
  
  
  res[[i]]  = list(nbmm_loglik   = modB2$fn(pars),
                   glmmTMB_loglik = modB2$fn(par2), 
                   params =  cbind(NBZIMM = pars, glmmTMB = par2),
                   AIC    =  2*(num_params - (-modB2$fn(pars)) ))
}

sum_a <- sum(sapply(res, function(x) x$AIC))


##############################################
##deseq

(deseq_mod$result$log2FoldChange)
ranef(rr_mod)






