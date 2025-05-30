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
param_index   =   1
source("simulation_analysis/initial_param0.R")
path = paste0("simulation_analysis/sim_data/",nsubj,"_",ntaxa,"/")
path
###########################################################
data       =   readRDS(paste0(path,"countdata.rds"))
metadata   =   readRDS(paste0(path,"metadata.rds"))
###########################################################
