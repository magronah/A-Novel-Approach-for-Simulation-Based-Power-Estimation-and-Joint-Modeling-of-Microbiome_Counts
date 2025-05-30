setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")

library(reformulas)
library(Matrix)
library(huge)
library(tidyverse)
library(dplyr)
library(MASS)
library(glmmTMB)
library(foreach)
library(NBZIMM)
library(nlme)
source("func2.R")
source("initial_param0.R")
####################################################################
ntaxa; nsubj; beta; betadisp; 

path  =  paste0("cov2/",nsubj,"_",ntaxa,"/")

if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE) 
  cat("Folder created at:", path, "\n")
} else {
  cat("Folder already exists at:", path, "\n")
}

####################################################################
dd     =   meta_data(ntaxa, nsubj)
dd$dummy  =    factor(1)
############Negative Binomial mixed models (NBMMs)############################
mod0      =     readRDS(paste0(path,"reference_model/nbmmref.rds"))

sim_list = list()

for(i in 1:ntaxa) {
  
  fit  =  mod0$fit[[i]]
  b    =  as.numeric(unlist(fit$coefficients$random))
  beta =  as.numeric(unlist((fit$coefficients$fixed)))
  theta_val  = as.numeric(VarCorr(fit)[, "StdDev"])[1]

  betadisp  =   as.numeric(unlist(fit$theta))
  sim_pars  =   lst(beta, betadisp, theta = log(theta_val), b)

  ddd       =   dd %>% filter(taxon == i)

  ####################################################################
  ## check if rr is being mapped correctly
  seed      =   i
  form      =   formula(paste0("taxon", i, "~ group + (1 | dummy)"))
  pars      =   simulate_new(RHSForm(form, as.form = TRUE), nsim = 1, seed = seed,
                             newdata    =  ddd,
                             newparams  =  sim_pars,
                             family     =  nbinom2,
                             return_val = "pars")
  
  stopifnot(b == pars[names(pars) == "b"])
  #####################################################################
  #now simulate data using b_est
  sim_list[[i]]  =  simulate_new(RHSForm(form, as.form = TRUE), nsim = nsim, seed = seed,
                      newdata   =  ddd,
                      newparams =  sim_pars,
                      family    =  nbinom2,
                      return_val = "sim")
  
}
names(sim_list)  =  paste0("taxon",1:ntaxa)

result <- lapply(1:nsim, function(j) {
  data.frame(do.call(cbind, lapply(sim_list, function(lst) lst[[j]])))
})

sim_dd  <- lapply(1:nsim, function(k) { 
  list(countdata =  result[[k]], 
       met_data   =  ddd)})

names(sim_dd)  =  paste0("sim", 1:nsim)
saveRDS(sim_dd, file=paste0(path,"sim_data/nbmm_otu_meta_list_withzi_taxa.rds"))
##############################################################################
# for(r in 1:ntaxa){
#   pp=sim_list[[r]][[1]]
#   stopifnot(sim_dd[[1]][[1]][,r] == pp)
# }

