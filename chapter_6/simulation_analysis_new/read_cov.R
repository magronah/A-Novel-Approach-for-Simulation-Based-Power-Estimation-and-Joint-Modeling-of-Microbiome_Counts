setwd("/home/agronahm/projects/def-bolker/agronahm/longitudinal_RR/")
library(foreach)
library(huge)
library(glmmTMB)
library(Matrix)
source("func.R")
############################################################
param_index  =  1
#########################################################
source("simulation_analysis_new/initial_param0.R")
path1  =  "simulation_analysis/sim_data/"

path = paste0("~/scratch/long/new_sim/",nsubj,"_",ntaxa,"_",ntime,"/coverage/rr")
path

files <-   list.files(path, full.names = TRUE)
res = foreach(i = files) %do%{
  readRDS(i)
}
length(res)
