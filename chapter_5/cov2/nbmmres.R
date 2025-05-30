setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")
library(MASS)
library(nlme)
library(NBZIMM)
library(foreach)
library(huge)
library(glmmTMB)
library(Matrix)
source("func2.R")
source("initial_param0.R")
############################################################
path = paste0("~/scratch/dataset/RR/coverage/",nsubj,"_",ntaxa,"/nbmm")
path


files =  list.files(path, full.names = TRUE)

res   = foreach(i = files,.packages = "NBZIMM") %do% {     
     mod   =   readRDS(i)
     pp    =   fixed(mod)$dist
     ppp   =  pp[(pp$variables) == "grouptreat",][["Estimate"]]
     names(ppp)  =  mod$response
     ppp
}

common_names <- Reduce(intersect, lapply(res, names))
filtered_res <- lapply(res, function(x) x[common_names])
res      <- do.call(cbind, filtered_res)

dd    =   as.data.frame(res)
#rownames(dd)  =  paste0("taxon",1:ntaxa)
colnames(dd)  =    paste0("nsim",1:ncol(dd))
saveRDS(dd, file = paste0("cov2/",nsubj,"_",ntaxa,"/nbmm.rds"))

