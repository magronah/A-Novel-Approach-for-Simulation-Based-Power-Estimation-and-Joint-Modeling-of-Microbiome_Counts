library(foreach)
library(huge)
library(glmmTMB)
library(Matrix)
library(here)

source("func2.R")
source("initial_param0.R")
############################################################
path = paste0("~/scratch/data/coverage/",nsubj,"_",ntaxa,"/confint/us")
path

files <-   list.files(path, full.names = TRUE)

res = foreach(i = files) %do% {
    readRDS(i)
}


names(res)  =  paste0("sim",1:length(files))
saveRDS(res, file = paste0(nsubj,"_",ntaxa,"/confint/us.rds"))




