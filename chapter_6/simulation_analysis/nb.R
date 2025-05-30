setwd("/home/agronahm/projects/def-bolker/agronahm/longitudinal_RR/")
library(RhpcBLASctl)
library(Matrix)
library(huge)
library(tidyverse)
library(dplyr)
library(NBZIMM)
library(DESeq2)
library(glmmTMB)
library(foreach)
library(doParallel)
###########################################################
source("func.R")
###########################################################
param_index   =   2
source("simulation_analysis/initial_param0.R")
path = paste0("simulation_analysis/sim_data/",nsubj,"_",ntaxa,"_",ntime,"/")
path
###########################################################
data       =   readRDS(paste0(path,"countdata.rds"))
metadata   =   readRDS(paste0(path,"metadata.rds"))
###########################################################
cc     =   commandArgs(trailingOnly  = TRUE)
i      =   as.integer(cc[1])

############################################################
ddd    =   data[[i]]
meta   =   metadata[[i]]
###########################################################
mod    =   mms(y = ddd, fixed = ~group*time  + offset(normalizer),
               random = ~ 1|subject,
               correlation = corAR1(form = ~ as.numeric(
                 factor(time, levels = unique(time)))| subject),
               data = meta,
               method = "nb",
               niter = 100)
###########################################################
file_path  =  paste0("~/scratch/long/",nsubj,"_",ntaxa,"_",ntime,"/nb/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}

saveRDS(mod, file=paste0(file_path,"mod",i,".rds"))




ddd    =   data[[i]]
meta   =   metadata[[i]]
###########################################################
mod    =   mms(y = ddd, fixed = ~group*time  + offset(normalizer),
                 random = ~ 1|subject,
                 correlation = corAR1(form = ~ as.numeric(
                   factor(time, levels = unique(time)))| subject),
                 zi_fixed = ~1,
                 data = meta,
                 method = "zinb",
                 niter = 100)
####################################################################
file_path  =  paste0("~/scratch/long/",nsubj,"_",ntaxa,"_",ntime,"/zinb/")
  
  if (!dir.exists(file_path)) {
    dir.create(file_path, recursive = TRUE)
    cat("Folder created at:", file_path, "\n")
  } else {
    cat("Folder already exists at:", file_path, "\n")
  }
  
saveRDS(mod, file=paste0(file_path,"mod",i,".rds"))
