library(glmmTMB)
library(ggplot2)
library(dplyr)
library(patchwork)
library(here)
##############################################################
atlass_path    <- paste0("real_data_analysis/atlass/results/")
pregnancy_path <- paste0("real_data_analysis/pregnancy/results/")
##############################################################
fig_path = "fig/"
source("func.R")
##############################################################
# Define filenames to load
mod_files <- c(
  "mod_rr.rds",
  "mod_us.rds")
#' confidence intervals
#' AIC calculation for each model
#' Run time
#' Statistical power
############################################################
# Load models 
atlass_models    <- load_models(atlass_path, mod_files)
pregnancy_models <- load_models(pregnancy_path, mod_files)

names(atlass_models)     =   c("RR", "US")
names(pregnancy_models)  =   c("RR", "US")
############################################################
path_vec  <-  c(atlass_path, pregnancy_path)

par_ctrl <- glmmTMBControl(
  parallel = list(n = 10, autopar = TRUE)
  # optCtrl  = list(eval.max=500, iter.max = 100)
)

gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")

#confidence intervals.
num  = 30
for(j in 1:2){
  j =1
  mod  <-  atlass_models[[j]]
  path <-  path_vec[j]
  for(i in 1:2){
    i =1
    name  =  names(atlass_models)[i]
    mmd  <-  mod
    df   <-  model.frame(mod) 
    dd    <-  df |> dplyr::rename(normalizer = "offset(normalizer)")
    lev  <-  leverage_brute(mmd, data = dd, inds = seq(num), eps = 0.1, 
                            fit_method = "update", scale = "response")
    
    saveRDS(lev, file= paste0(path,"leverage", name, ".rds"))
  }
  
  
}



# 