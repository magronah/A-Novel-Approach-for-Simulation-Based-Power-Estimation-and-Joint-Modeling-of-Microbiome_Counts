setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")
library(glmmTMB)
library(ggplot2)
library(dplyr)
library(Matrix)
##############################################################
path1   =   paste0("real_data/")
source(paste0(path1,"fun.R"))
#########################################################
filenames <- c("rr_mod.rds")
autism_path <- paste0("real_data/autism_data/results/")
autism_models <- load_models(autism_path, filenames)
#########################################################
mod   =  autism_models[[1]] 
#########################################################
par_ctrl <- glmmTMBControl(parallel = list(n = 10, autopar = TRUE),
                           optCtrl  = list(eval.max=100, iter.max = 10)
                           )

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
file_path  =  paste0(autism_path,"leverage/rr/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}

saveRDS(lev, file= paste0(file_path,"lev",i, ".rds"))
