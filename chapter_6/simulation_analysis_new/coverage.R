setwd("/home/agronahm/projects/def-bolker/agronahm/longitudinal_RR/")
library(RhpcBLASctl)
library(Matrix)
library(huge)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(glmmTMB)
library(foreach)
###########################################################
source("func.R")
###########################################################
param_index   =   2
source("simulation_analysis_new/initial_param0.R")
path = paste0("~/scratch/long/new_sim/",nsubj,"_",ntaxa,"_",ntime,"/us")
path0 = paste0("simulation_analysis_new/sim_data/",nsubj,"_",ntaxa,"_",ntime,"/")
path

grp_label = "grouptreatment:time"
pp= zinbmm_confint(mod, mean_count = mean_list[[1]],
                   group_label = grp_label,
                   conf_level = .95)

leftjoin(pp, true_param, by = "param_name")

###########################################################
data       =   readRDS(paste0(path0,"countdata.rds"))
true_param =   readRDS(paste0(path0,"results/true_param.rds"))

mean_list  =   lapply(data,function(x){colMeans(x)})

files  =   list.files(path, full.names = TRUE)
#cc     =   commandArgs(trailingOnly  = TRUE)
#i      =   as.integer(cc[1])
#i   =  1
###########################################################
res = foreach(i = 1:length(files)) %do% {
mean_count    =   mean_list[[i]]
mod           =   readRDS(files[i])
confint       =   wald_confint3(mod, ntaxa =  ntaxa,
                                mean_count =  mean_count,
                                true_param =  true_param,                                     
                                mod_name   =  NULL)
confint

}

saveRDS(res,file = paste0(path0,"results/coverage_us.rds"))




if(FALSE){
confint       =   wald_confint3(mod, ntaxa =  ntaxa,
                                mean_count =  mean_count,
                                true_param =  true_param,
                                mod_name   =  NULL)
###########################################################
file_path  =  paste0("~/scratch/long/new_sim/",nsubj,"_",ntaxa,"_",ntime,"/coverage/us/")
if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}
saveRDS(confint,file = paste0(file_path,"confint_",i,".rds"))
}
