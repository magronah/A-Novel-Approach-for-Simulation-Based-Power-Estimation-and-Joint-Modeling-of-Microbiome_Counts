setwd("/home/agronahm/projects/def-bolker/agronahm/longitudinal_RR/")
library(RhpcBLASctl)
library(Matrix)
library(huge)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(glmmTMB)
library(foreach)
###########################################################
source("func.R")
#######################################
param_index   =   2
source("simulation_analysis_new/initial_param0.R")
path = paste0("simulation_analysis_new/sim_data/",nsubj,"_",ntaxa,"_",ntime,"/")
path
###########################################################
data   =   readRDS(paste0(path,"dd_long.rds"))
cc     =   commandArgs(trailingOnly  = TRUE)
i      =   as.integer(cc[1])


dat   =   data[[i]]
###########################################################
form2	   =   count ~ 1  +
        (group*time|taxon) +
       (1|nugget) +
      (1 | subject:taxon) 

# try with   different
  par_ctrl <- glmmTMBControl(
    parallel = list(n = 10, autopar = TRUE)
  )

  gprior  <- data.frame(prior = "gamma(2, 2.5)",
                        class = "theta_sd",
                        coef = "")

  options(glmmTMB_openmp_debug = TRUE)
  blas_set_num_threads(1)

  system.time(
    fit  <-  glmmTMB(form2, data  = dat,
                     family  = poisson,
                     ziformula  = ~1,
                     prior   = gprior,
                     REML    = TRUE,
                     control = par_ctrl
    )
  )

file_path  =  paste0("~/scratch/long/new_sim/",nsubj,"_",ntaxa,"_",ntime,"/us/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}

saveRDS(fit, file=paste0(file_path,"mod",i,".rds"))
