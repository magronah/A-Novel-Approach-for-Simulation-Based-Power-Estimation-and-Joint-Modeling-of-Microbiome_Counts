setwd("/home/agronahm/projects/def-bolker/agronahm/longitudinal_RR/")
library(glmmTMB)
library(ggplot2)
library(dplyr)
library(patchwork)
library(nlme)
library(MASS)
library(Matrix)
library(NBZIMM)
library(huge)
library(foreach)
library(dplyr)
#############################################
source("func.R")
################################################
cc     =   commandArgs(trailingOnly  = TRUE)
j      =   as.integer(cc[1])
#####################################################
param_index   =   j
source("simulation_analysis_new/initial_param0.R")
path = paste0("~/scratch/long/new_sim/",nsubj,"_",ntaxa,"_",ntime,"/zinb")
path0 = paste0("simulation_analysis_new/sim_data/",nsubj,"_",ntaxa,"_",ntime,"/")
path
###########################################################
data	   =   readRDS(paste0(path0,"countdata.rds"))
true_param =   readRDS(paste0(path0,"results/true_param.rds"))

mean_list  =   lapply(data,function(x){colMeans(x)})

files  =   list.files(path, full.names = TRUE)

grp_label = "grouptreatment:time"

res = foreach(i = 1:length(files), .errorhandling = "remove") %do% {
mod           =    readRDS(files[i])
mean_mu       =    mean_list[[i]]
true_para     =    true_param[true_param$param_name  %in% mod$responses,]
mean_count    =    mean_mu[names(mean_list[[i]]) %in% mod$responses]

confint       =    zinbmm_confint(mod, mean_count = mean_count,
                              group_label = grp_label,
                              conf_level = .95)

left_join(confint, true_param, by = "param_name")

}

saveRDS(res,file = paste0(path0,"results/coverage_zinb.rds"))
