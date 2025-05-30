setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(foreach)
source("load_glmmTMB.R")

path =   "~/scratch/dataset/new_sim/coverage/rr_mod/mod"

cc   =   commandArgs(trailingOnly  = TRUE)
n    =   as.integer(cc[1]) 

step_size =  200 

start  =   (n - 1) * step_size + 1
end    =   (n * step_size)  

res = foreach(i =  start:end, .combine = "cbind") %do% {
  mod  =   readRDS(paste0(path,i,".rds"))
  dd   =  (ranef(mod)$cond$taxon$grouptreat)
  names(dd)  =  paste0("taxon",1:length(dd))
  dd
}

saveRDS(res, file = paste0("~/scratch/dataset/new_sim/coverage/rr_res/param",n,".rds"))



#for(i in start:end){
#  mod  =   readRDS(paste0(path,i,".rds"))
#  dd   =  (ranef(mod)$cond$taxon$grouptreat)
#names(dd)  =  paste0("taxon",1:length(dd))
#saveRDS(dd, file = paste0("~/scratch/dataset/new_sim/coverage/rr_res/res",i, ".rds"))
#}


