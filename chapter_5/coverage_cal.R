library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(scam)
library(ggrastr)
library(metR)
library(latex2exp)
source("reproducible/new_sim/func.R")
########################################################################################
path        =   paste0(getwd(), "/reproducible/new_sim/data/")
fig_path    =   paste0(getwd(), "/reproducible/new_sim/fig/")
###################################################################
mean_count  =   readRDS(paste0(path,"/mean_count_dd.rds"))
true_param  =  readRDS(paste0(path, "true_param_new.rds"))
rr_est      =  readRDS(paste0(path, "rr_parametric_new.rds"))
us_est      =  readRDS(paste0(path, "us_only_parametric.rds"))
deseq_est   =  readRDS(paste0(path, "deseq_parametric2_new.rds"))
nbmm_est   =  readRDS(paste0(path, "nbmm_parametric.rds"))
##################################################################
lmean_abund   =   rowMeans(mean_count)
lfoldchange   =   true_param$true_param
names(lmean_abund)  =   rownames(lfoldchange)
########################################################################## 
##convert to long format
rr_dd       =   dd_long(rr_est,true_param,label="with_rr")
us_dd       =   dd_long(us_est,true_param,label="without_rr")
deseq_dd    =   dd_long(deseq_est,true_param,label="deseq")
nbmm_dd     =   dd_long(nbmm_est,true_param,label="nbmm")
########################################################################### 
rr_para_conf      =   para_confint(rr_dd, true_param, alpha = 0.05)
us_para_conf      =   para_confint(us_dd, true_param, alpha = 0.05)
deseq_para_conf   =   para_confint(deseq_dd, true_param, alpha = 0.05)
nbmm_para_conf   =    para_confint(nbmm_dd, true_param, alpha = 0.05)
################################################################################
pval_rr    = coverage_cal(rr_para_conf)$coverage
pval_us    = coverage_cal(us_para_conf)$coverage
pval_deseq = coverage_cal(deseq_para_conf)$coverage
pval_nbmm  = coverage_cal(nbmm_para_conf)$coverage
##################################################################
### fit GAM
fit_rr     =   gam_fit2(pval_rr,lfoldchange,lmean_abund)
fit_us     =   gam_fit2(pval_us,lfoldchange,lmean_abund)
fit_deseq  =   gam_fit2(pval_deseq,lfoldchange,lmean_abund)
fit_nbmm   =   gam_fit2(pval_nbmm,lfoldchange,lmean_abund)
###############################################################
combined_dd_rr          =   fit_rr$combined_data
power_estimate_rr       =   fit_rr$power_estimate

combined_dd_us         =   fit_us$combined_data
power_estimate_us      =   fit_us$power_estimate

combined_dd_deseq      =   fit_deseq$combined_data
power_estimate_deseq   =   fit_deseq$power_estimate

combined_dd_nbmm       =   fit_nbmm$combined_data
power_estimate_nbmm    =   fit_nbmm$power_estimate
############################################################
cont_breaks  = c(0.8,0.975) #c(0.8,0.9,0.95,0.975) # seq(0,1,0.01)
pp_rr     =  contour_plot_fun(combined_dd_rr,power_estimate_rr,cont_breaks)
pp_us     =  contour_plot_fun(combined_dd_us,power_estimate_us,cont_breaks)
pp_deseq  =  contour_plot_fun(combined_dd_deseq,power_estimate_deseq,cont_breaks)
pp_nbmm   =  contour_plot_fun(combined_dd_nbmm,power_estimate_nbmm,cont_breaks)

(pp_rr|pp_deseq|pp_nbmm) +  plot_layout(guides = "collect")


