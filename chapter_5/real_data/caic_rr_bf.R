setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")
library(ggplot2)
library(dplyr)
library(patchwork)
library(AICcmodavg)
source("reproducible/leverage_funs.R")
library(glmmTMB)
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
  "rrzi_each_mod.rds"
  )
######################################################################
# Load models for autism data
autism_path <- paste0("real_data/autism_data/results/")
soil_path   <- paste0("real_data/soil_data/results/")

# Load models for autism data
autism_models <- load_models(autism_path, filenames)
soil_models   <- load_models(soil_path, filenames)

# Assigning names 
mod_names <- c("RR","RRzi")
names(autism_models)    =   names(soil_models)  =    mod_names
mod_list    =   lst(autism_models,soil_models)
rr_list <- lapply(mod_list, function(x){x$RR})
nam   =  c("autism_data","soil_data")
#####################################################
cc =   commandArgs(trailingOnly  = TRUE)
i         =   as.integer(cc[1])
i =1
mod   =  rr_list[[i]]
dd    =  model.frame(mod)
n <- nrow(dd)
batch_size <- 50  
grp_list <- split(1:n, ceiling(seq_along(1:n) / batch_size))
length(grp_list)

model =  fm1; H_matrice = list()
for(j in 1:length(grp_list)){
  H_matrix <- leverage_brute_force(fm1, data = dd,grp_index = j,
                                         grp_list = grp_list)
  sum(diag(H_matrix))
}

condlik3 <- sum(dnbinom(model.response(model.frame(mod)),
                        mu = fitted(mod), size = sigma(mod),
                        log = TRUE))
#####################################################
res =  c(clik = condlik3)

file_path  =  paste0("real_data/",nam[i],"/results/")
if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}

saveRDS(lever, file=paste0(file_path,"leverage_rrzi.rds"))
saveRDS(res, file=paste0(file_path,"caic_rrzi.rds"))

