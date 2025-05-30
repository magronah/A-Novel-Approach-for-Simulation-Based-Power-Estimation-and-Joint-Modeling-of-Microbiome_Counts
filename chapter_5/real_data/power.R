library(tidyverse)
library(foreach)
library(doParallel)
library(here)
#################################################
path1   =   paste0("real_data/")
source(paste0(path1,"/fun.R"))
#################################################
# Load models for autism data
autism_path <- paste0("real_data/autism_data/results/")
atlass_path <- paste0("real_data/atlass_data/results/")
crohn_path <- paste0("real_data/CrohnD_data/results/")
soil_path   <- paste0("real_data/soil_data/results/")
#########################################################
# Define filenames to load
filenames <- c(
  "rr_mod.rds",
  "rrzi_each_mod.rds",
  "us_mod.rds",
  "uszi_each_mod.rds",
  "nbmm_aicc.rds",
  "zinbmm_aicc.rds",
  "deseq_aicc.rds"  
)
#################################################
# Load models for autism data
autism_models <- load_models(autism_path, filenames)
atlass_models <- load_models(atlass_path, filenames)
crohn_models <- load_models(crohn_path, filenames)
soil_models   <- load_models(soil_path, filenames)

# Assigning names 
mod_names <- c("RR","RRzi","US","USzi","Nbmm","Zinbmm","DE")
names(autism_models)    =   names(atlass_models)  =    mod_names
names(crohn_models)    =   names(soil_models)    =    mod_names
mod_list    =   lst(autism_models,atlass_models,
                    crohn_models, soil_models)
#
pr = autism_models$US
(pr$sdr$diag.cov.random)
