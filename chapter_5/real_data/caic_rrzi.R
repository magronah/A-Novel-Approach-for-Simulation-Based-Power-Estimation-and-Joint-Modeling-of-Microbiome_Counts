setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")
#library(glmmTMB)
library(ggplot2)
library(dplyr)
library(patchwork)
library(AICcmodavg)
source("reproducible/leverage_funs.R")
library(glmmTMB, lib.loc = "/home/agronahm/projects/def-bolker/agronahm/glmmTMB_lev")
library(RTMB)
library(Matrix)
library(lme4)
library(cAIC4)
##############################################################
path1   =   paste0("real_data/")
source(paste0(path1,"fun.R"))
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

######################################################################
# Load models for autism data
autism_path <- paste0("real_data/autism_data/results/")
atlass_path <- paste0("real_data/atlass_data/results/")
crohn_path  <- paste0("real_data/CrohnD_data/results/")
soil_path   <- paste0("real_data/soil_data/results/")

# Load models for autism data
autism_models <- load_models(autism_path, filenames)
atlass_models <- load_models(atlass_path, filenames)
crohn_models <- load_models(crohn_path, filenames)
soil_models   <- load_models(soil_path, filenames)

# Assigning names 
mod_names <- c("RR","RRzi","US","USzi","NB","ZNB","DE")
names(autism_models)    =   names(atlass_models)  =    mod_names
names(crohn_models)    =   names(soil_models)    =    mod_names
mod_list    =   lst(autism_models,atlass_models,crohn_models,soil_models)
rr_list <- lapply(mod_list, function(x){x$RRzi})
nam   =  c("autism_data","atlass_data","CrohnD_data","soil_data")
#####################################################
cc =   commandArgs(trailingOnly  = TRUE)
i         =   as.integer(cc[1])

mod   =  rr_list[[i]]
condlik3 <- sum(dnbinom(model.response(model.frame(mod)),
                        mu = fitted(mod), size = sigma(mod),
                        log = TRUE))

system.time(lever <- leverage(mod))
trace_hat <-   sum(lever)
cdf3	  <-   trace_hat+1
#####################################################
res =  c(clik = condlik3, cdf = cdf3, caic = 2*(-condlik3 + cdf3))

file_path  =  paste0("real_data/",nam[i],"/results/")
if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}

saveRDS(lever, file=paste0(file_path,"leverage_rrzi.rds"))
saveRDS(res, file=paste0(file_path,"caic_rrzi.rds"))

