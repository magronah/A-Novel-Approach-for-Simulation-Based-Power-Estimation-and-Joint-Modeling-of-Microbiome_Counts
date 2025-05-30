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
mean_count_rr  =   readRDS(paste0(path,"sim_data/rr_mean_count.rds"))
res_rr	       =   readRDS(paste0(path,"parametric_res/rr.rds"))
##########################################################
lfc_rr       =   rowMeans(res_rr)
#lmean       =   log(mean)
pval_rr      =   pvalue_cal(res_rr)

saveRDS(pval_rr, file  =  paste0(path,"pvalues/rr.rds"))
###########################################################
fit_rr       =    gam_fit(pval_rr, lfc_rr, mean_count_rr,
                          grid_len = 100, alpha_level = 0.05)

saveRDS(fit_rr, file  =  paste0(path,"gam_fit/rr.rds"))

#######################rrzi###############################
mean_count_rrzi   =    readRDS(paste0(path,"sim_data/rrzi_mean_count.rds"))
res_rrzi          =    readRDS(paste0(path,"parametric_res/rrzi.rds"))
##########################################################
lfc_rrzi	  =    rowMeans(res_rrzi)
#lmean     =   log(mean)
pval_rrzi	  =    pvalue_cal(res_rrzi)

saveRDS(pval_rrzi, file  =  paste0(path,"pvalues/rrzi.rds"))
###########################################################
fit_rrzi      =    gam_fit(pval_rrzi,lfc_rrzi, mean_count_rrzi,
                            grid_len = 100, alpha_level = 0.05)

saveRDS(fit_rrzi, file  =  paste0(path,"gam_fit/rrzi.rds"))


