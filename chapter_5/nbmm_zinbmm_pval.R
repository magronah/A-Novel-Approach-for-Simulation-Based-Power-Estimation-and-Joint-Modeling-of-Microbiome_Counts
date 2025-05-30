setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")

library(foreach)
library(doParallel)
library(NBZIMM)

library(dplyr)



##############################################################
### nbmm
directory_path <-   "~/scratch/dataset/new_sim/nbmm_pval/"
files          <-   list.files(directory_path, full.names = TRUE)


names_list   =  foreach(i = files) %do% {
      dd     =    readRDS(i); names(dd)
}

common_names <- Reduce(intersect, names_list)

print(names_list)

res   = foreach(i = files, .combine="cbind") %do% {
    tryCatch({
    pval_dd   =  readRDS(i)
    pval      =  pval_dd[names(pval_dd) %in% common_names]
    stopifnot(names(pval) == common_names)
    pval
   }, error = function(e) {
      message(paste("Error in file:", i))
      return(NA)  # Return NA or another value that makes sense in your context
    })
  }

print(res)
res            =    as.data.frame(res)
rownames(res)  =    common_names
colnames(res)  =    paste0("nsim",1:ncol(res))
saveRDS(res, file = paste0(getwd(), "/reproducible/new_sim/data/pval_nbmm.rds"))


##############################################################
if(FALSE){
### zinbmm
directory_path2 <-   "~/scratch/dataset/new_sim/zinbmm_pval/"
files          <-   list.files(directory_path2, full.names = TRUE)


names_list   =  foreach(i = files) %do% {
      dd     =    readRDS(i); names(dd)
}

common_names <- Reduce(intersect, names_list)


res   = foreach(i = files, .combine="cbind") %do% {
    tryCatch({
    pval_dd   =  readRDS(i)
    pval      =  pval_dd[names(pval_dd) %in% common_names]
    stopifnot(names(pval) == common_names)
    pval
   }, error = function(e) {
      message(paste("Error in file:", i))
      return(NA)  # Return NA or another value that makes sense in your context
    })
  }

res            =    as.data.frame(res)
rownames(res)  =    common_names
colnames(res)  =    paste0("nsim",1:ncol(res))
saveRDS(res, file = paste0(getwd(), "/reproducible/new_sim/data/pval_zinbmm.rds"))
}
