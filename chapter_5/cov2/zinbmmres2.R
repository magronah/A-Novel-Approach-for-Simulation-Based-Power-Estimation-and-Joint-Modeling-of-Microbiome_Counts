setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")

library(foreach)
library(huge)
library(glmmTMB)
library(Matrix)
source("func2.R")
source("initial_param0.R")
############################################################
path = paste0("~/scratch/dataset/RR/coverage/",nsubj,"_",ntaxa,"/zinbmm2")
path

files <-   list.files(path, full.names = TRUE)

res = foreach(i = files, .combine ="cbind") %do% {
  mod  =   readRDS(i)
  lapply(mod, function(x){summary(x)$coefficients$cond["grouptreat", "Estimate"]})
}


dd            =  as.data.frame(res)
rownames(dd)  =  paste0("taxon",1:nrow(dd))
colnames(dd)  =  paste0("sim",1:ncol(dd))
saveRDS(dd, file = paste0("cov2/",nsubj,"_",ntaxa,"/zinbmm.rds"))


