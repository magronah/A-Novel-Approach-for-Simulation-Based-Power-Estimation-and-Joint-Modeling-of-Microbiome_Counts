library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
source("reproducible/new_sim/func.R")
path        =   paste0(getwd(), "/reproducible/new_sim/")
#########################################################
data_10_50   =  load_data("10_50")
#data_50_100  =  load_data("50_100")
data_100_300 =  load_data("100_300")
########################################################
errors_fun <- function(est_list) {
  
  true_param  =  est_list$true_param
  nbmm_est    =  est_list$nbmm_est
  zinbmm_est  =  est_list$zinbmm_est
  deseq_est   =  est_list$deseq_est
  us_est      =  est_list$us_est
  rr_est      =  est_list$rr_est
  
  err = list(nbmm_error    =    error_cal(nbmm_est, true_param, model = "nbmm"),
             zinbmm_error  =    error_cal(zinbmm_est, true_param, model = "zinbmm"),
             deseq_error   =    error_cal(deseq_est, true_param, model = "deseq"),
             us_error      =    error_cal(us_est, true_param, model = "without_rr"))
  
  avg_mse    =  do.call(rbind, lapply(err, function(x){x$avg_mse}))
  avg_bias   =  do.call(rbind, lapply(err, function(x){x$avg_bias}))
  rownames(avg_mse) =  rownames(avg_bias) =  NULL
  
  n             =    nrow(avg_bias) 
  rr_error      =    error_cal(rr_est, true_param, model = "rr")
  mse_rr        =    rr_error$avg_mse
  bias_rr       =    rr_error$avg_bias
  
  mse_rr     =   do.call(rbind, replicate(n, mse_rr, simplify = FALSE))
  bias_rr    =   do.call(rbind, replicate(n, bias_rr, simplify = FALSE))
  
  avg_mse$type  =    mse_rr$type = paste0("mse", 1:nrow(avg_mse))
  avg_bias$type =    bias_rr$type = paste0("bias", 1:nrow(avg_bias))
  
  avg_mse       =    rbind(avg_mse,  mse_rr)
  avg_bias      =    rbind(avg_bias, bias_rr)
  
  list(avg_mse = avg_mse, avg_bias = avg_bias)
}

pp       =   errors_fun(data_10_50)
mse_dd   =   pp$avg_mse
bias_dd  =   pp$avg_bias

n=15
ppm = ggplot(mse_dd, aes(x =  type, y = average_value,group = type,  
                         color = factor(model, levels = c("rr", "without_rr","nbmm","deseq","zinbmm")))) +
  geom_point(size=2) +
  scale_color_manual(values = c("rr" = "blue", "without_rr" = "red", "deseq" = "green",
                                "nbmm" = "orange", "zinbmm" = "purple")) +
  labs(color = "model") +
  ylab("average mse across taxa") +
  custom_theme(n) +
  theme(axis.title.x = element_blank())

#ggsave(paste0(fig_path,"avgmse.png"), plot = ppm, width = 8, height = 6, dpi = 500)
################################################################################
ppb = ggplot(bias_dd, aes(x =  type, y = average_value,group = type,  
                          color = factor(model, levels = c("rr", "without_rr","deseq","zinbmm","nbmm")))) +
  geom_point(size=2) +
  scale_color_manual(values = c("rr" = "blue", "without_rr" = "red", "deseq" = "green",
                                "nbmm" = "orange", "zinbmm" = "purple")) +
  labs(color = "model") +
  ylab("average bias across taxa") +
  custom_theme(n) +
  theme(axis.title.x = element_blank())

#ggsave(paste0(fig_path,"avgbias.png"), plot = ppb, width = 8, height = 6, dpi = 500)
