setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(foreach)
source("load_glmmTMB.R")

path =   "~/scratch/dataset/new_sim/us_only/mod"

cc   =   commandArgs(trailingOnly  = TRUE)
i    =   as.integer(cc[1])

mod  =   readRDS(paste0(path,i,".rds"))


dd   =  (ranef(mod)$cond$taxon$grouptreat)

names(dd)  =  paste0("taxon",1:length(dd))
saveRDS(dd, file = paste0("~/scratch/dataset/new_sim/us_only_res/res",i, ".rds"))






