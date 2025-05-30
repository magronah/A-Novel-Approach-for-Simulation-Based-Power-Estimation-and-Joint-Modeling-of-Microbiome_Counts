setwd("/home/agronahm/projects/def-bolker/agronahm/longitudinal_RR/")
library(glmmTMB)
library(ggplot2)
library(dplyr)
library(patchwork)
library(here)
##############################################################
path <- paste0("real_data_analysis/pregnancy/results/")
##############################################################
fig_path = "fig/"
source("func.R")
##############################################################
mod  <-  readRDS(paste0(path, "mod_us.rds"))

par_ctrl <- glmmTMBControl(parallel = list(n = 10, autopar = TRUE))

gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")
############################################################
cc     =   commandArgs(trailingOnly  = TRUE)
i      =   as.integer(cc[1])
############################################################
mmd  <-  mod
df   <-  model.frame(mod) 
dd   <-  df |> dplyr::rename(normalizer = "offset(normalizer)")

lev  <-  leverage_brute_modified(mmd, data = dd, inds = i, eps = 0.1,
                        fit_method = "update", scale = "response",
                        pred_method  = "predict",
                        progress = TRUE)
############################################################
file_path  =  paste0(path,"leverage/us/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}

saveRDS(lev, file= paste0(file_path,"lev",i, ".rds"))

