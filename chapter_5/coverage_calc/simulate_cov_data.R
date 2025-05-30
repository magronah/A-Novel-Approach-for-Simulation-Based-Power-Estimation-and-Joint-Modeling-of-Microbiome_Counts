setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
source("reproducible/new_sim/Load_Packages.R")
source("reproducible/new_sim/func.R")
source("reproducible/new_sim/param.R")
source("reproducible/new_sim/initial_param0.R")
source("load_glmmTMB.R")
path = paste0(getwd(),"/reproducible/new_sim/coverage_calc/data/")
###################################################################
ntaxa=100; nsubj=50; beta; betadisp; theta_true; betazi; nsim    =  20000
form <- count ~ 1 + us(1 + group|taxon) + us(0 + taxon | subject)
#####################################################################
dd          =   met_data(ntaxa, nsubj)
true_pars0  =   lst(beta, theta = theta_true, betadisp)
####################################################################  
#### simulate fixed bs
pars0 <- simulate_new(RHSForm(form, as.form = TRUE), nsim = 1, seed = seed,
                     newdata   =  dd,
                     newparams =  true_pars0,
                     family     =  nbinom2, 
                     return_val = "pars"
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
                         seed      =  seed)

res_dd1 = foreach(i = 1:length(count_list))  %do% {
  dd$count = count_list[[i]]
  dd
}
names(res_dd1)  =   paste0("sim", 1:nsim)

####################################################################  
#res_dd1   =  readRDS(paste0(path,"sim_count_list_nozi.rds"))
#### convert to otu table
res_dd2   =  otu_meta_lst_fun(res_dd1)
####################################################################  
saveRDS(true_b, file = paste0(path,"true_param2.rds"))


saveRDS(res_dd1, file = paste0(path,"sim_count_list_nozi2.rds"))
saveRDS(res_dd2, file = paste0(path,"otu_meta_list_nozi.rds"))
