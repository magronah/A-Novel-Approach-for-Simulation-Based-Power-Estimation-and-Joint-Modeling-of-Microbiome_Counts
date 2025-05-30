library(RhpcBLASctl)
library(Matrix)
library(huge)
library(glmmTMB)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(here)
###########################################################
source("func2.R")
source("initial_param0.R")
###########################################################
path = paste0("~/scratch/data/",nsubj,"_",ntaxa,"/rr/")
path
###########################################################
cc    =   commandArgs(trailingOnly  = TRUE)
i     =   as.integer(cc[1])
###########################################################
data	=   readRDS(paste0(path,"mod",i,".rds"))
###########################################################
path  =   paste0("~/scratch/data/coverage/",nsubj,"_",ntaxa,"/rr")
