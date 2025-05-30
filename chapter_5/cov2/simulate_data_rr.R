setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")

library(reformulas)
library(Matrix)
library(huge)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(glmmTMB)
library(foreach)
source("func2.R")
source("initial_param0.R")
###################################################################
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
form   =   count ~ 1 + us(1 + group|taxon) +  rr(0 + taxon | subject,2)
####################################################################
## simulate from rr and rrzi models
mod   =  list(rr     =    readRDS(paste0(path,"rrref.rds")),
              rrzi   =    readRDS(paste0(path,"rrziref.rds")))

 
theta_true1 =  c(
  get_theta_logSD(n_gt, seed = seed),
  get_theta_corr(n_gt, n_gt,  seed=seed),
  get_theta_rr(d = 2, n = ntaxa, logsdvec = log(c(2, 1)))
)


for(i in 1:length(mod)){
      b     =   mod[[i]]$obj$report()$b
  sim_pars  =   lst(b, beta, theta = theta_true1, betadisp)
  
  #######################################################################
  ## check if rr is being mapped correctly
  pars      =   simulate_new(RHSForm(form, as.form = TRUE), nsim = 1, seed = seed,
                             newdata    =  dd,
                             newparams  =  sim_pars,
                             family     =  nbinom2,
                             return_val = "pars")
  
  stopifnot(b == pars[names(pars) == "b"])
  stopifnot(theta_true1 == pars[names(pars)=="theta"])
  ################################################################################
  #now simulate data using b_est
  count_list <- simulate_new(RHSForm(form, as.form = TRUE),
                             newdata   = dd,
                             newparams = sim_pars,
                             family = nbinom2,
                             nsim  = nsim,
                             seed = seed)
  
  res_dd1 = foreach(j = 1:length(count_list))  %do% {
    dd$count = count_list[[j]]
    dd
  }
  names(res_dd1)  =   paste0("sim", 1:nsim)
  
  ####################################################################  
  #### convert to otu table
  res_dd2   =  otu_meta_lst_fun(res_dd1)
  ####################################################################  
  nam  =  names(mod)[[i]]
  saveRDS(res_dd1, file = paste0(path,"sim_data/", nam,"_sim_count_list_nozi.rds"))
  saveRDS(res_dd2, file = paste0(path,"sim_data/", nam, "_otu_meta_list_nozi.rds"))
  ####################################################################  
  
   sim_pars1          =   c(sim_pars, thetazi =  log(diff(qlogis(c(0.12, 0.92)))/4))
  sim_pars1$betazi    =   qlogis(((0.12) +(0.92))/2)
  
  #########################################################################
  #### Sanity check
  pars3 <- simulate_new(RHSForm(form, as.form = TRUE), nsim = 1, seed = NULL,
                        newdata    =  dd,
                        newparams  =  sim_pars1,
                        family     =  nbinom2, 
                        ziformula  = ~1 + (1|taxon),
                        return_val = "pars"
  )
  
  stopifnot(theta_true1 ==  pars3[names(pars3)=="theta"])
  stopifnot(b ==  pars3[names(pars3)=="b"])
  
  ################## simulate count data zero inflation per taxa
  ztaxa_count_list <- simulate_new(RHSForm(form, as.form = TRUE),
                                   newdata   =  dd,
                                   newparams =  sim_pars1,
                                   family    =  nbinom2,
                                   ziformula = ~1 + (1|taxon),
                                   nsim      =  nsim,
                                   seed      =  seed)
  
  ztaxa_res_dd1 = foreach(i = 1:length(ztaxa_count_list))  %do% {
    dd$count = ztaxa_count_list[[i]]
    dd
  }
  
  names(ztaxa_res_dd1)  =   paste0("sim", 1:nsim)
  ####################################################################  
  #### convert to otu table
  ztaxa_res_dd2   =  otu_meta_lst_fun(ztaxa_res_dd1)
  ####################################################################  
  saveRDS(ztaxa_res_dd1, file = paste0(path,"sim_data/",nam,"_sim_count_list_withzi_taxa.rds"))
  saveRDS(ztaxa_res_dd2, file = paste0(path,"sim_data/",nam,"_otu_meta_list_withzi_taxa.rds"))
  
  mean_count  =   colMeans(do.call(rbind,lapply(ztaxa_res_dd2, function(x) (x$countdata))))
  saveRDS(mean_count, file = paste0(path,"sim_data/",nam,"_mean_count.rds"))
}



