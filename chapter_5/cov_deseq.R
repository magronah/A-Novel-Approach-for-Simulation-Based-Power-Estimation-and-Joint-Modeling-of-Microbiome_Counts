library(foreach)
library(huge)
library(glmmTMB)
library(Matrix)
library(here)

source("func2.R")
source("initial_param0.R")
############################################################
path = paste0("~/scratch/data/",nsubj,"_",ntaxa,"/deseq")
path

files <-   list.files(path, full.names = TRUE)

res = foreach(i = files) %do% {
    mod =  readRDS(i)
    deseq_wald_confint(mod)
}


names(res)  =  paste0("sim",1:length(files))
saveRDS(res, file = paste0(nsubj,"_",ntaxa,"/confint/deseq.rds"))




