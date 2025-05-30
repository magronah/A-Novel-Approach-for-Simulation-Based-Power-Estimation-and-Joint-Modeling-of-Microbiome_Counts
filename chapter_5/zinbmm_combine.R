setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")

library(foreach)
library(doParallel)
library(NBZIMM)


directory_path <-   "~/scratch/dataset/zinbmm_mod/"

cc    =   commandArgs(trailingOnly  = TRUE)
i     =   as.integer(cc[1])

mod   =   readRDS(paste0(directory_path,"mod",i,".rds"))

nsim  =   length(mod$fit)
print(length(mod$fit))

res   = foreach(j = 1:nsim, .combine="c", .packages = "NBZIMM") %do% {
    tryCatch({
         coef(mod$fit[[j]])$fixed[["grouptreat"]]
    }, error = function(e) {
      message(paste("Error in file:", i, "at simulation:", j, "-", e))
      return(NA)
    })
  }


names(res)  =  names(mod$fit)
saveRDS(res, file = paste0("~/scratch/dataset/zinbmm_res/res",i,".rds" ))


