library(foreach)
library(huge)
library(NBZIMM)
library(nlme)
library(glmmTMB)
library(Matrix)
library(here)

source("func2.R")
source("initial_param0.R")
############################################################
path = paste0("~/scratch/data/",nsubj,"_",ntaxa,"/nbmm")
path

files <-   list.files(path, full.names = TRUE)

res = foreach(i = files) %do% {
    mod =  readRDS(i)
    zinbmm_confint(mod)
}


names(res)  =  paste0("sim",1:length(files))
saveRDS(res, file = paste0(nsubj,"_",ntaxa,"/confint/nbmm.rds"))



