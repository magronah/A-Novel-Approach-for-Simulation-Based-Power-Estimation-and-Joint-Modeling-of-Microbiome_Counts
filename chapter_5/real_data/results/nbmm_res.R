setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")

library(foreach)
library(doParallel)
library(NBZIMM)

path  =   "~/scratch/dataset/new_sim/nbmm_mod/"
mod   =   readRDS(paste0(path,"mod.rds"))


cc    =   commandArgs(trailingOnly  = TRUE)


nsim  =   length(mod$fit) 
res   =   foreach(j = 1:nsim, .combine="c", .packages = "NBZIMM") %do% {
   tryCatch({
     coef(mod$fit[[j]])[["grouptreat"]]
   }, error = function(e) {
     message(paste("Error in file:", i, "at simulation:", j, "-", e))
     return(NA)
   })
 }
 
names(res)  =  names(mod$fit)
saveRDS(res, file = paste0("~/scratch/dataset/new_sim/autis_nbmm_res.rds" ))
 

