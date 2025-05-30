library(huge)
library(glmmTMB)
library(Matrix)
library(here)
###########################################################
source("func2.R")
source("initial_param0.R")
###########################################################
mm    =  "rrzi"
path  =  paste0("~/scratch/data/",nsubj,"_",ntaxa,"/",mm)
path
files  <-   list.files(path, full.names = TRUE)
############################################################
path_dir  =   paste0(nsubj,"_",ntaxa,"/")
cc        =   commandArgs(trailingOnly  = TRUE)
i         =   as.integer(cc[1])

###########################################################
model     =   readRDS(files[i])
mu_count  =   readRDS(paste0(path_dir,"mean_count.rds"))
nam	  =   paste0("mod",i)
###########################################################
path_d     =  paste0("~/scratch/data/coverage/",nsubj,"_",ntaxa)

save_path  =  paste0(path_d,"/sdreport/",mm,"/")


if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
  cat("Folder created at:", save_path, "\n")
} else {
  cat("Folder already exists at:", save_path, "\n")
}



sdr_path  =  paste0(path_d,"/fullsdr/",mm,"/")


if (!dir.exists(sdr_path)) {
  dir.create(sdr_path, recursive = TRUE)
  cat("Folder created at:", sdr_path, "\n")
} else {
  cat("Folder already exists at:", sdr_path, "\n")
}

mod_nam  = paste0("mod",i)
confint   <-  wald_confint(mod   =  model,
                           ntaxa =  ntaxa,
                           mean_count =  mu_count,
                           mod_name =  mod_nam,
                           sdreport_path =  save_path, 
                           fullsdr_path  =  sdr_path)


confint_path  =  paste0(path_d,"/confint/",mm,"/")

if (!dir.exists(confint_path)) {
  dir.create(confint_path, recursive = TRUE)
  cat("Folder created at:", confint_path, "\n")
} else {
  cat("Folder already exists at:", confint_path, "\n")
}

saveRDS(confint, file = paste0(confint_path,"/mod",i,".rds"))


