setwd("/home/agronahm/projects/def-bolker/agronahm/longitudinal_RR/")
library(foreach)
library(huge)
library(glmmTMB)
library(Matrix)
source("func.R")
############################################################
param_index  =  2
#########################################################
source("simulation_analysis_new/initial_param0.R")
path1  =  "simulation_analysis_new/sim_data/"

path = paste0("~/scratch/long/new_sim/",nsubj,"_",ntaxa,"_",ntime,"/us")
path

files <-   list.files(path, full.names = TRUE)

res = foreach(i = files, .combine ="cbind") %do% {
    mod  =   readRDS(i)
    dd   =   ranef(mod,condVar = FALSE)$cond$taxon$`grouptreatment:time`
    dd
}


dd            =  as.data.frame(res)
rownames(dd)  =  paste0("taxon",1:nrow(dd))

sim_numbers <- gsub("\\D+", "", basename(files))
colnames(dd)  =  paste0("sim", sim_numbers)
saveRDS(dd, file = paste0(path1,nsubj,"_",ntaxa,"_",ntime,"/results/us.rds"))

