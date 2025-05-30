setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")
library(scam)
library(RhpcBLASctl)
library(Matrix)
library(huge)
library(glmmTMB)
library(tidyverse)
library(dplyr)
library(DESeq2)
source("func2.R")
source("initial_param0.R")
path = paste0("cov2/",nsubj,"_",ntaxa,"/")
path
###########################################################
mean_count_us  =   readRDS(paste0(path,"sim_data/us_mean_count.rds"))
res_us         =   readRDS(paste0(path,"parametric_res/us.rds"))
##########################################################
lfc_us       =   rowMeans(res_us)
#lmean       =   log(mean)
pval_us      =   pvalue_cal(res_us)

saveRDS(pval_us, file  =  paste0(path,"pvalues/us.rds"))
###########################################################
fit_us       =    gam_fit(pval_us, lfc_us, mean_count_us,
                          grid_len = 100, alpha_level = 0.05)

saveRDS(fit_us, file  =  paste0(path,"gam_fit/us.rds"))

#######################rrzi###############################
mean_count_uszi   =    readRDS(paste0(path,"sim_data/uszi_mean_count.rds"))
res_uszi          =    readRDS(paste0(path,"parametric_res/uszi.rds"))
##########################################################
lfc_uszi          =    rowMeans(res_uszi)
#lmean     =   log(mean)
pval_uszi         =    pvalue_cal(res_uszi)

saveRDS(pval_uszi, file  =  paste0(path,"pvalues/uszi.rds"))
###########################################################
fit_uszi      =    gam_fit(pval_uszi,lfc_uszi, mean_count_uszi,
                            grid_len = 100, alpha_level = 0.05)

saveRDS(fit_uszi, file  =  paste0(path,"gam_fit/uszi.rds"))

