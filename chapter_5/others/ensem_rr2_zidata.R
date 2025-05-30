setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
source("reproducible/new_sim/func.R")
source("reproducible/new_sim/param.R")
source("reproducible/new_sim/Load_Packages.R")
source("reproducible/new_sim/initial_param0.R")
source("load_glmmTMB.R")
path = paste0(getwd(),"/reproducible/new_sim/data/")
###########################################################
data      =   readRDS(paste0(path,"sim_count_list_withzi.rds"))
form      =   count ~ 1 + us(1 + group|taxon) +  rr(0 + taxon | subject,2)
################################################################
cc      =   commandArgs(trailingOnly  = TRUE)
i       =   as.integer(cc[1]) 
dd      =   data[[i]]
################################################################
##now add normalization constant 
df  =  otu_meta_fun(dd)

##Fit the model
par_ctrl <- glmmTMBControl(
  parallel = list(n = 10, autopar = TRUE)
)

gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")

form2  =   update(form,.~. + offset(normalizer))

system.time(
  fit  <-  glmmTMB(form2, data = df,
                   family  = nbinom2, 
                   prior   = gprior,
                   REML    = TRUE,
                   control = par_ctrl
  )
)

saveRDS(fit, file=paste0("~/scratch/dataset/new_sim/rr_mod2zi/mod", i, ".rds"))

