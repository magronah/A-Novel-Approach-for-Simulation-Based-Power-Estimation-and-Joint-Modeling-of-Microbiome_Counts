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
path = paste0("cov2/",nsubj,"_",ntaxa,"/")
path
###########################################################
data      =   readRDS(paste0(path,"sim_data/","rr_sim_count_list_withzi_taxa.rds"))
form      =   count ~ 1 + us(1 + group|taxon) +  rr(0 + taxon | subject,2)
################################################################
cc =   commandArgs(trailingOnly  = TRUE)
i         =   as.integer(cc[1])

##Fit the model
par_ctrl <- glmmTMBControl(
  parallel = list(n = 10, autopar = TRUE)
)


 dd     =   data[[i]]
  ################################################################
  ##now add normalization constant 
  df  =  otu_meta_fun(dd)

  gprior  <- data.frame(prior = "gamma(2, 2.5)",
                        class = "theta_sd",
                        coef = "")

  form2  =   update(form,.~. + offset(normalizer))

  options(glmmTMB_openmp_debug = TRUE)
  blas_set_num_threads(1)

mod  <- tryCatch({
       glmmTMB(form2, data = df,
                     family  = nbinom2,
                     prior   = gprior,
                     REML    = TRUE,
                     control = par_ctrl)
}, error = function(e) {
     message("Error in first attempt, trying again without prior...")
       glmmTMB(form2, data = df,
                     family  = nbinom2,
                     REML    = TRUE,
                     control = par_ctrl)
})


file_path  =  paste0("~/scratch/dataset/RR/coverage","/",nsubj,"_",ntaxa,"/","rr/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}

saveRDS(mod, file=paste0(file_path,"mod",i,".rds")) 

