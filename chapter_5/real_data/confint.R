library(tidyverse)
library(foreach)
library(doParallel)
library(Matrix)
library(here)
#################################################
path1   =   paste0("real_data/")
source(paste0(path1,"fun.R"))
#################################################
# Load models for autism data
autism_path <- paste0("real_data/autism_data/results/")
atlass_path <- paste0("real_data/atlass_data/results/")
crohn_path <- paste0("real_data/CrohnD_data/results/")
soil_path   <- paste0("real_data/soil_data/results/")
#########################################################
# Define filenames to load
filenames <- c(
  "rr_mod.rds",
  "rrzi_each_mod.rds",
  "us_mod.rds",
  "uszi_each_mod.rds",
  "nbmm_mod.rds",
  "zinbmm_mod.rds",
  "deseq_mod_est.rds"  
)
#################################################
# Load models for autism data
autism_models <- load_models(autism_path, filenames)
atlass_models <- load_models(atlass_path, filenames)
crohn_models <- load_models(crohn_path, filenames)
soil_models   <- load_models(soil_path, filenames)

# Assigning names 
mod_names <- c("RR","RRzi","US","USzi","NB","ZNB","DE")
names(autism_models)    =   names(atlass_models)  =    mod_names
names(crohn_models)    =   names(soil_models)    =    mod_names
mod_list    =   lst(autism_models,atlass_models,
                    crohn_models, soil_models)
###########################################################
## us and rr models only
us_rr_models <- lapply(mod_list, function(x){
  x[grep("^(RR|US)", names(x))]
})
###########################################################
autism_us_rr  =   us_rr_models$autism_models
atlass_us_rr  =   us_rr_models$atlass_models
crohn_us_rr   =   us_rr_models$crohn_models
 soil_us_rr   =   us_rr_models$soil_models
 ############################################################
ntaxa_file    <-  "ntaxa.rds"
autism_ntaxa  <-  load_models(autism_path, ntaxa_file)
atlass_ntaxa  <-  load_models(atlass_path, ntaxa_file)
crohn_ntaxa   <-  load_models(crohn_path,  ntaxa_file)
soil_ntaxa    <-  load_models(soil_path, ntaxa_file) 
############################################################
datasets <- list(
  autism  = list(model = autism_us_rr, ntaxa = autism_ntaxa, path  = autism_path),
  atlass  = list(model = atlass_us_rr, ntaxa = atlass_ntaxa, path = atlass_path),
  crohn   = list(model = crohn_us_rr,  ntaxa = crohn_ntaxa, path = crohn_path),
  soil    = list(model = soil_us_rr,   ntaxa = soil_ntaxa, path = soil_path)
)

group_vec =  c(autism = "groupNT", atlass = "groupobese",
               crohn =  "groupCD", 
               soil = "groupT")
#name = "soil"

for (name in names(datasets)) {
   model <- datasets[[name]]$model
   ntax <- as.numeric(datasets[[name]]$ntaxa)
   path_dir <- datasets[[name]]$path
  
   mu_count  =   readRDS(paste0(path_dir,"mean_count.rds"))
   grp_lab   =   group_vec[[name]]
   confint_us <- wald_confint3(mod  = model$US,
                              ntaxa = ntax,
                              mean_count = mu_count,
                              mod_name = "us_wald3",
                              path = path_dir,
                              grp_label = grp_lab)

   confint_uszi <- wald_confint3(mod  = model$USzi,
                              ntaxa = ntax,
                              mean_count = mu_count,
                              mod_name = "uszi_wald3",
                              path = path_dir,
                              grp_label = grp_lab)

   saveRDS(confint_us, file = paste0(path_dir, "CI_us_wald3.rds"))
   saveRDS(confint_uszi, file = paste0(path_dir, "CI_uszi_wald3.rds"))
   ###########################################################
   confint_rr <- wald_confint3(mod  = model$RR,
                              ntaxa = ntax,
                              mean_count = mu_count,
                              mod_name = "rr_wald3",
                              path = path_dir,
                              grp_label = grp_lab)

   saveRDS(confint_rr, file = paste0(path_dir, "CI_rr_wald3.rds"))
   # ###########################################################
   confint_rrzi <- wald_confint3(mod  = model$RRzi,
                                ntaxa = ntax,
                                mean_count = mu_count,
                                mod_name = "rrzi_wald3",
                                path = path_dir,
                                grp_label = grp_lab)

   saveRDS(confint_rrzi, file = paste0(path_dir, "CI_rrzi_wald3.rds"))
}

#' what was the problem with fitting rrzi? 
#' simulation studies coverage and confidence interval and power
###########################################################
## deseq model and rr models only
datasets_deseq <- list(
  autism  = list(model = autism_models, path  = autism_path),
  atlass  = list(model = atlass_models, path = atlass_path),
  crohn   = list(model = crohn_models, path = crohn_path),
  soil    = list(model = soil_models, path = soil_path)
)

for(i in 1:length(datasets_deseq)){
  deseq_est  =   datasets_deseq[[i]]$model$DE
  path      <-   datasets_deseq[[i]]$path
  mu_count  =   readRDS(paste0(path,"mean_count.rds"))
  
  CI_deseq   =   deseq_wald_confint(deseq_est, mean_count = mu_count)
  saveRDS(CI_deseq, file = paste0(path, "CI_deseq.rds"))
}
######################################################
i = 1
for(i in 1:length(datasets_deseq)){
  mod_nbmm   <-   datasets_deseq[[i]]$model$Nbmm
  mod_zinbmm <-   datasets_deseq[[i]]$model$Zinbmm
  #################################################
  path       <-   datasets_deseq[[i]]$path
  mu_count   <-   readRDS(paste0(path,"mean_count.rds"))
  #################################################
  name_vec   <-  mod_nbmm$variables$dist
  grp_label  <-  name_vec[grep("group", name_vec)]
  #################################################
  CI_nbmm   <- zinbmm_confint(mod_nbmm, mean_count = mu_count,
                            group_label = grp_label,
                            conf_level = .95)
  
 plot(CI_nbmm$est_param)
 plot(CI_zinbmm$est_param)
 
 CI_zinbmm  <- zinbmm_confint(mod_zinbmm, mean_count = mu_count,
                              group_label = grp_label,
                              conf_level = .95)
  
 #################################################  
 saveRDS(CI_nbmm, file = paste0(path, "CI_nbmm.rds"))
 saveRDS(CI_zinbmm, file = paste0(path, "CI_zinbmm.rds"))
}
