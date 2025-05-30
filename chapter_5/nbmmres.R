library(NBZIMM)
library(foreach)
library(huge)
library(glmmTMB)
library(Matrix)
source("func2.R")
source("initial_param0.R")
############################################################
path = paste0("~/scratch/dataset/RR","/",nsubj,"_",ntaxa,"/nbmm")
path


files =   list.files(path, full.names = TRUE)

res   = foreach(i = files,.combine = "cbind",.packages = "NBZIMM") %do% {     
     mod   =   readRDS(i)
     pp    =   fixed(mod)$dist
     pp[(pp$variables) == "grouptreat",][["Estimate"]]
}
dd    =   as.data.frame(res)
rownames(dd)  =  paste0("taxon",1:ntaxa)
colnames(dd)  =    paste0("nsim",1:ncol(dd))
saveRDS(dd, file = paste0(getwd(),"/",nsubj,"_",ntaxa,"/nbmm.rds"))

