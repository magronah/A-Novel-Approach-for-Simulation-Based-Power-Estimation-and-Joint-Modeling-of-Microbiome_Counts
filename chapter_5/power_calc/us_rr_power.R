library(tibble)
library(DESeq2)
library(Matrix)
library(dplyr)
library(huge)
library(foreach)
library(scam)
library(tidyverse)
#############################################
source("func2.R")
source("initial_param0.R")

#source("/project/6006158/agronahm/reduced_rank_mm/func2.R")
#source("/project/6006158/agronahm/reduced_rank_mm/initial_param0.R")

# Define sets of nsubj and ntaxa
parameter_sets <- list(
 list(nsubj = 100, ntaxa = 300),
 list(nsubj = 150, ntaxa = 500),
 list(nsubj = 200, ntaxa = 600)
)


#models = list()
# Loop through parameter sets
for (params in parameter_sets) {
  nsubj <- params$nsubj
  ntaxa <- params$ntaxa
  
  
  # Set file paths
  # path <- paste0("/project/6006158/agronahm/reduced_rank_mm/", nsubj, "_", ntaxa, "/")
  path <- paste0(nsubj, "_", ntaxa, "/")
  
  # Load data and extract parameters
  dd <- load_data(path)
  effect_size <- dd$dd$true_param$true_param
  pval_data <- dd$dd[c("rr", "rrzi", "us", "uszi")]
 
  # Compute p-values for all models
  pvals <- lapply(pval_data, pvalue_cal)
  names(pvals) <- c("rr","rrzi", "us", "uszi")
  
  # Read mean count data
  mean_count <- readRDS(paste0(path, "mean_count.rds"))

  # Fit GAM models for each p-value set
  models <- lapply(pvals, gam_fit, effect_size, mean_count, grid_len = 500, alpha_level = 0.05)
  names(models) <- c("rr","rrzi", "us", "uszi")
  
  # # Create output directory if it doesn't exist
  file_path <- paste0(path, "GAM/")
  if (!dir.exists(file_path)) {
    dir.create(file_path, recursive = TRUE)
    cat("Folder created at:", file_path, "\n")
  } else {
    cat("Folder already exists at:", file_path, "\n")
  }
  
  # Save models
  lapply(names(models), function(name) {
     saveRDS(models[[name]], file = paste0(file_path, name, ".rds"))
   })
  # 
   # Create output directory if it doesn't exist
   file_path2 <- paste0(path, "pvalues/")
   if (!dir.exists(file_path2)) {
     dir.create(file_path2, recursive = TRUE)
     cat("Folder created at:", file_path2, "\n")
   } else {
     cat("Folder already exists at:", file_path2, "\n")
   }
   
  # # Save pvalues
   lapply(names(pvals), function(name) {
     saveRDS(pvals[[name]], file = paste0(file_path2, name, ".rds"))
   })
}

