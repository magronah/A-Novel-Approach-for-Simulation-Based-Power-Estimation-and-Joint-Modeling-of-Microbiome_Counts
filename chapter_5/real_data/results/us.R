#setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
path        =   paste0(getwd(),"/reproducible/new_sim/real_data/")
##############################################################
## have nodes =  1 and 
# cpu-per-task= more
library(tidyr)
library(tidyverse)
library(tibble)
library(DESeq2)
library(RhpcBLASctl)
source("load_glmmTMB.R")
#library(glmmTMB)
source(paste0(path,"fun.R"))
source(paste0(path,"autism/load_data.R"))
##############################################################
dd_long     =   df_long(dd_filt)
dd_long     =   left_join(dd_long, met_dd, by ="subject")  
###########################################################
dds        =   DESeqDataSetFromMatrix(dd_filt,met_dd, ~group)
dds        =   DESeq(dds,sfType ="poscounts",minReplicatesForReplace=Inf) 
normalizer =   data.frame(normalizer = sizeFactors(dds)) %>% 
  rownames_to_column("subject")

df      =    left_join(dd_long, normalizer, by ="subject")  
###########################################################
form	  =   count ~ 1 + us(1 + group|taxon) + offset(normalizer)

##Fit the model
par_ctrl <- glmmTMBControl(
  parallel = list(n = 6, autopar = TRUE),
  optArgs  = list(maxit = 10000)
)

gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")

options(glmmTMB_openmp_debug = TRUE)
blas_set_num_threads(1)

tt=system.time(
  fit  <-  glmmTMB(form, data = df,
                   family  = nbinom2, 
                   prior   = gprior,
                   REML    = TRUE,
                   control = par_ctrl
  )
)

saveRDS(fit, file=paste0(path,"autism/data/us_mod.rds"))
saveRDS(tt, file=paste0(path,"autism/data/us_runtime",".rds"))



