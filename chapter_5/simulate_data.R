library(reformulas)
library(RhpcBLASctl)
library(Matrix)
library(huge)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(glmmTMB)
library(foreach)
library(here)
source("func2.R")
source("initial_param0.R")
###################################################################
ntaxa; nsubj; beta; betadisp; 

path  =  paste0(nsubj,"_",ntaxa,"/")

if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE) 
  cat("Folder created at:", path, "\n")
} else {
  cat("Folder already exists at:", path, "\n")
}

## specify the model
form <- count ~ 1 + us(1 + group|taxon) + us(0 + taxon | subject)
#####################################################################
dd          =   meta_data(ntaxa, nsubj)
true_pars0  =   lst(beta, theta = theta_true, betadisp)

p=get_theta_corr(ntaxa, nsubj,  seed=seed)
plot(p)
####################################################################  
#### simulate fixed bs
pars0 <- simulate_new(RHSForm(form, as.form = TRUE), nsim = 1, seed = seed,
                     newdata    =  dd,
                     newparams  =  true_pars0,
                     family     =  nbinom2, 
                     return_val =  "pars"
)

stopifnot(theta_true == pars0[names(pars0)=="theta"])
b_true_all      =   pars0[names(pars0)=="b"]
true_b_index    =   seq(2,(2*ntaxa),2)
true_param      =   b_true_all[true_b_index]
true_b          =   data.frame(param_name = paste0("taxon", 1:ntaxa), 
                           true_param = true_param)
################################################################################
#now simulate with bs fixed to bval
true_pars1  =  c(true_pars0, lst(b = b_true_all))
pars1 <- simulate_new(RHSForm(form, as.form = TRUE), nsim = 1, seed = NULL,
                      newdata    =  dd,
                      newparams  =  true_pars1,
                      family     =  nbinom2, 
                      return_val = "pars"
)
stopifnot(b_true_all  == pars1[names(pars1) == "b"])
####################################################################  
#### simulate count data no zero inflation
count_list <- simulate_new(RHSForm(form, as.form = TRUE),
                         newdata   =  dd,
                         newparams =  true_pars1,
                         family    =  nbinom2,
                         nsim      =  nsim,
                         seed   =  seed)

res_dd1 = foreach(i = 1:length(count_list))  %do% {
  dd$count = count_list[[i]]
  dd
}
names(res_dd1)  =   paste0("sim", 1:nsim)
####################################################################  
#### convert to otu table
res_dd2   =  otu_meta_lst_fun(res_dd1)
####################################################################  
saveRDS(true_b, file = paste0(path,"true_param.rds"))
saveRDS(res_dd1, file = paste0(path,"sim_count_list_nozi.rds"))
saveRDS(res_dd2, file = paste0(path,"otu_meta_list_nozi.rds"))
####################################################################  
true_pars3          =   c(true_pars1, thetazi =  log(diff(qlogis(c(0.12, 0.92)))/4))
set.seed(seed)
true_pars3$betazi   =  qlogis(((0.12) +(0.92))/2)
#########################################################################
#### Sanity check
pars3 <- simulate_new(RHSForm(form, as.form = TRUE), nsim = 1, seed = NULL,
                     newdata    =  dd,
                     newparams  =  true_pars3,
                     family     =  nbinom2, 
                     ziformula  = ~1 + (1|taxon),
                     return_val = "pars"
)
stopifnot(theta_true ==  pars3[names(pars3)=="theta"])
stopifnot(b_true_all ==  pars3[names(pars3)=="b"])

################## simulate count data zero inflation per taxa
ztaxa_count_list <- simulate_new(RHSForm(form, as.form = TRUE),
                            newdata   = dd,
                            newparams = true_pars3,
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
saveRDS(ztaxa_res_dd1, file = paste0(path,"sim_count_list_withzi_taxa.rds"))
saveRDS(ztaxa_res_dd2, file = paste0(path,"otu_meta_list_withzi_taxa.rds"))

mean_count  =   colMeans(do.call(rbind,lapply(ztaxa_res_dd2, function(x) (x$countdata))))
saveRDS(mean_count, file = paste0(path, "mean_count.rds"))
