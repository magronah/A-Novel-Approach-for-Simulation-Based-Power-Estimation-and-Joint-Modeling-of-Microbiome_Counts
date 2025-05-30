setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(foreach)
#source("load_glmmTMB.R")
devtools::load_all("~/Documents/PhD/glmmTMB/glmmTMB")

library(RhpcBLASctl)

path =   "~/sharcnet_mount/new_sim/100_300/rr_mod/mod"
#readRDS(path)

#path =   "~/scratch/dataset/new_sim/100_300/rr_mod/mod"

#cc   =   commandArgs(trailingOnly  = TRUE)
#i    =   as.integer(cc[1]) 

i    =   1
mod  =   readRDS(paste0(path,i,".rds"))

system.time(
  blas_set_num_threads(1),
  ranef(mod)
)

coef(mod)


options(glmmTMB_openmp_debug = TRUE)
dd   =  (ranef(mod)$cond$taxon$grouptreat)

names(dd)  =  paste0("taxon",1:length(dd))
saveRDS(dd, file = paste0("~/scratch/dataset/new_sim/100_300/rr_mod_res/res",i,".rds"))

# saveRDS(dd, file = paste0("~/scratch/dataset/new_sim/100_300/rr_mod_res/res",i,".rds"))



