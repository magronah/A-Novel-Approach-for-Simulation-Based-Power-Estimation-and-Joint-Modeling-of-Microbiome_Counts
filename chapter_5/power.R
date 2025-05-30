setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
source("reproducible/new_sim/func.R")
###############################################################################
path        =   paste0(getwd(), "/reproducible/new_sim/data/")
fig_path    =   paste0(getwd(), "/reproducible/new_sim/fig/")
###############################################################################
true_param  =   readRDS(paste0(path, "true_param_new.rds"))
deseq_est   =   readRDS(paste0(path, "deseq_parametric_new.rds"))
###############################################################################
deseq_est   =  readRDS(paste0(path, "deseq_parametric_new.rds"))
View(deseq_est)
