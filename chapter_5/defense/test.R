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
source("defense/initial_param0.R")
###################################################################
path  =  paste0("defense/",nsubj,"_",ntaxa,"/")

if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE) 
  cat("Folder created at:", path, "\n")
} else {
  cat("Folder already exists at:", path, "\n")
}
###################################################################
form <- count ~ 1 + us(1 + group|taxon) + us(0 + taxon | subject)
#####################################################################
dd          =   meta_data(ntaxa, nsubj)
s =  seq(1,5,length=nsim)

ztaxa_count_list = list()
for(j in 1:nsim){
  ppp  =    get_theta_corr1(ntaxa, nsubj,  seed=seed)

  theta_true =  c(
    get_theta_logSD(n_gt, seed = seed),
    get_theta_corr1(n_gt, n_gt,  seed=seed),
    get_theta_logSD(ntaxa, seed = seed),
    ppp*s[[j]]
  )
print(range(ppp*s[[j]]))
true_pars0  =   lst(beta, theta = theta_true, betadisp)
###################################################################
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

ztaxa_count_list[[j]] <- simulate_new(RHSForm(form, as.form = TRUE),
                                 newdata   = dd,
                                 newparams = true_pars3,
                                 family    =  nbinom2,
                                 ziformula = ~1 + (1|taxon),
                                 nsim      =  1,
                                 seed      =  seed)[[1]]
}


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


