setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")
library(RhpcBLASctl)
library(Matrix)
library(huge)
library(glmmTMB)
library(tidyverse)
library(dplyr)
library(DESeq2)
source("func2.R")
source("initial_param0.R")
path  =  paste0("cov2/",nsubj,"_",ntaxa,"/")
path
####################################################################
data	  =   readRDS(paste0(path,"sim_data/","zinbmm_otu_meta_list_withzi_taxa.rds"))
################################################################
cc	=   commandArgs(trailingOnly  = TRUE)
i	  =   as.integer(cc[1])

par_ctrl <- glmmTMBControl(
  parallel = list(n = 10, autopar = TRUE)
)

dd	=   data[[i]]
################################################################
countdata  =   dd$countdata
met_dd     =   dd$met_data

met_dd$dummy = factor(1)

otu_count  =   t(countdata)
dds        =   DESeqDataSetFromMatrix(otu_count,met_dd, ~group)
dds        =   DESeq(dds,sfType ="poscounts",minReplicatesForReplace=Inf)
normalizer =   sizeFactors(dds)


options(glmmTMB_openmp_debug = TRUE)
blas_set_num_threads(1)


gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")

#lapply(mod, function(x){summary(x)$coefficients$cond["grouptreat", "Estimate"]})

mod_list  =   list()

for(j in 1:ntaxa){
  
  met_dd$count   =   countdata[,j]
  
  fit  <- tryCatch({
    mod_fit  =  glmmTMB(count ~ group + (1 | dummy) + offset(normalizer), 
                  data = met_dd,
                  family  = nbinom2,
                  ziformula  = ~1,
                  prior   = gprior,
                  REML    = TRUE,
                  control = par_ctrl)
  }, error =  function(e){
    message("Error in first attempt, trying again without prior...")
    mod_fit  =  glmmTMB(count ~ group + (1 | dummy) + offset(normalizer), 
                  data = met_dd,
                  family  = nbinom2,
                  ziformula  = ~1,
                  REML    = TRUE,
                  control = par_ctrl)
    
  })
  
  mod_list[[j]]  = mod_fit
  
}

file_path  =  paste0("~/scratch/dataset/RR/coverage","/",nsubj,"_",ntaxa,"/","zinbmm2/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}

names(mod_list)  =  colnames(countdata)
saveRDS(mod_list, file=paste0(file_path,"mod",i,".rds"))



