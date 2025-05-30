library(foreach)
library(huge)
library(glmmTMB)
library(Matrix)
source("func2.R")
source("initial_param0.R")

path = paste0("~/scratch/dataset/RR","/",nsubj,"_",ntaxa,"/rrzi")
path
##########################################################
files <-   list.files(path, full.names = TRUE)

res = foreach(i = files, .combine ="cbind") %do% {
      tryCatch({
    mod  =   readRDS(i)
    dd   =  (ranef(mod,condVar = FALSE)$cond$taxon$grouptreat)
    dd}, error = function(e) {
    message("Error in file: ", i, " - ", e)
    NULL
  })
}

dd            =  as.data.frame(res)
rownames(dd)  =  paste0("taxon",1:nrow(dd))
colnames(dd)  =  paste0("sim",1:ncol(dd))

saveRDS(dd, file = paste0(getwd(),"/",nsubj,"_",ntaxa,"/rrzi.rds"))
