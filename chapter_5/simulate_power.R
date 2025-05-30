setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
source("reproducible/simp_single/func.R")
source("reproducible/simp_single/param.R")
source("reproducible/simp_single/Load_Packages.R")
source("reproducible/simp_single/initial_param0.R")
source("load_glmmTMB.R")
path = paste0(getwd(),"/reproducible/simp_single/data/")
##################################################################
form  <- count ~ 1 + us(1 + group|taxon) + us(0 + taxon | subject)
meta_dd     =   met_data(ntaxa, nsubj)
true_pars   =   lst(beta, theta = theta_true, betadisp, betazi)
#####################################################################
## simulate  count data and effect
nsim_       =  50     
effects_dd  =  sim_effect(form, meta_dd, true_pars, nsim=nsim_)
count_dd    =  sim_count(effects_dd,theta_true,meta_dd)
##################################################################
res = foreach(i  =   1:length(count_dd)) %do% {
  count_ff      =   count_dd[[i]]
  otu_tab_meta  =    spread(count_ff, key = taxon, value = count) 
  colnames(otu_tab_meta)  =  c("subject","group", paste0("taxon",1:ntaxa))
  meta_dd     =   otu_tab_meta %>% dplyr::select(c(subject, group))
  
otu_table     =   otu_tab_meta %>% 
                    dplyr::select(-c(subject, group)) %>%
                    `rownames<-`(meta_dd$subject)

 lst(otu_table,meta_dd)
}

names(res) = paste0("sim",1:length(count_dd))

saveRDS(effects_dd, file = paste0(path,"effect_size_dd.rds"))
saveRDS(count_dd,   file = paste0(path,"count_dd_pow.rds"))
saveRDS(res,   file = paste0(path,"otu_meta_dd_pow.rds"))
### save mean abundance here
